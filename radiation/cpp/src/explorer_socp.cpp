/*
 * Copyright (c) 2015, The Regents of the University of California (Regents).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the name of the copyright holder nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Please contact the author(s) of this library if you have any questions.
 * Author: David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 */

///////////////////////////////////////////////////////////////////////////////
//
// Exploration on a 2D grid. Tries to find the specified number of radiation
// sources (located at random lattice points) by choosing trajectories of
// the specified number of steps that maximize mutual information between
// simulated measurements and the true map.
//
///////////////////////////////////////////////////////////////////////////////

#include <explorer_socp.h>

#ifdef SYSTEM_OSX
#include <gurobi_c++.h>
#endif

#include <glog/logging.h>
#include <random>
#include <string>
#include <math.h>

namespace radiation {

// Constructor/destructor.
ExplorerSOCP::~ExplorerSOCP() {}
ExplorerSOCP::ExplorerSOCP(unsigned int num_rows, unsigned int num_cols,
                           unsigned int num_sources, double regularizer,
                           unsigned int num_steps, double fov,
                           unsigned int num_samples,
                           double epsilon,
                           const std::vector<Source2D>& sources,
                           const GridPose2D& initial_pose)
  : Explorer2D(num_rows, num_cols, num_sources, regularizer, fov,
               sources, initial_pose),
    num_steps_(num_steps),
    num_samples_(num_samples),
    epsilon_(epsilon) {}
ExplorerSOCP::ExplorerSOCP(unsigned int num_rows, unsigned int num_cols,
                           unsigned int num_sources, double regularizer,
                           unsigned int num_steps, double fov,
                           unsigned int num_samples, double epsilon)
  : Explorer2D(num_rows, num_cols, num_sources, regularizer, fov),
    num_steps_(num_steps),
    num_samples_(num_samples),
    epsilon_(epsilon) {}

// Plan a new trajectory.
bool ExplorerSOCP::PlanAhead(std::vector<GridPose2D>& trajectory) {
#ifdef SYSTEM_LINUX
  VLOG(1) << "Your machine does not have Gurobi. Sorry.";
  return false;
#endif

#ifdef SYSTEM_OSX
  // Use GUROBI to set up and solve a SOCP. Since GUROBI uses
  // exceptions to communicate errors, we enclose the entire planner in
  // a try/catch block.
  try {
    // Create an empty environment and model.
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    // Generate conditional entropy vector.
    Eigen::VectorXd hzx;
    std::vector<unsigned int> trajectory_ids;
    map_.GenerateEntropyVector(num_samples_, num_steps_, pose_, fov_,
                               hzx, trajectory_ids);
    CHECK(hzx.rows() == trajectory_ids.size());

    // Create GUROBI variables.
    std::vector<GRBVar> variables;
    for (unsigned int ii = 0; ii < trajectory_ids.size(); ii++) {
      GRBVar p = model.addVar(0.0, /* lower bound */
                              1.0, /* upper bound */
                              0.0, /* objective coefficient */
                              GRB_CONTINUOUS, /* type */
                              "p" + std::to_string(ii) /* name */);
      variables.push_back(p);
    }

    GRBVar slack = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "slack");

    // Set up objective function, probability constraint, and slack constraint.
    GRBLinExpr objective = 0.0;
    GRBLinExpr prob_constraint = 0.0;
    GRBQuadExpr slack_constraint = 0.0;
    for (unsigned int ii = 0; ii < variables.size(); ii++) {
      objective += hzx(ii) * variables[ii];
      prob_constraint += variables[ii];
      slack_constraint += variables[ii] * variables[ii];
    }

    // Add slack to the objective.
    objective -= epsilon_ * slack;

    // Add these to the model.
    model.setObjective(objective, GRB_MAXIMIZE);
    model.addConstr(prob_constraint == 1.0, "probability constraint");
    model.addQConstr(slack_constraint <= slack * slack, "slack constraint");

    // Optimize the model.
    model.optimize();

    // Sample from the optimal distribution at random. Choose a random value in
    // [0.0, 1.0] and step through the distribution until we exceed this value.
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    const double random_value = unif(rng);
    double cdf = 0.0;
    unsigned int vertex_id = 0;
    while (cdf < random_value && vertex_id < variables.size())
      cdf += variables[vertex_id++].get(GRB_DoubleAttr_X);

    const unsigned int trajectory_id = trajectory_ids[--vertex_id];

    // Decode this trajectory id.
    trajectory.clear();
    DecodeTrajectory(trajectory_id, num_steps_, pose_, trajectory);
    return true;
  } catch (GRBException exception) {
    VLOG(1) << "Gurobi error code : " << exception.getErrorCode();
    VLOG(1) << exception.getMessage();
  } catch (...) {
    VLOG(1) << "Exception during optimization.";
  }

  return false;
#endif
}

} // namespace radiation
