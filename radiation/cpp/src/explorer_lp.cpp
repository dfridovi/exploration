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

#include <explorer_lp.h>

#include <GLUT/glut.h>
#include <glog/logging.h>
#include <gurobi_c++.h>
#include <random>
#include <string>
#include <math.h>

namespace radiation {

// Constructor/destructor.
ExplorerLP::~ExplorerLP() {}
ExplorerLP::ExplorerLP(unsigned int num_rows, unsigned int num_cols,
                       unsigned int num_sources, double regularizer,
                       unsigned int num_steps, double fov,
                       unsigned int num_samples)
  : map_(num_rows, num_cols, num_sources, regularizer),
    num_steps_(num_steps),
    num_samples_(num_samples),
    pose_(0.0, 0.0, 0.0),
    fov_(fov) {
  // Set up a random number generator.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, num_rows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, num_cols - 1);
  std::uniform_real_distribution<double> unif_angle(0.0, 2.0 * M_PI);

  // Choose random sources.
  for (unsigned int ii = 0; ii < num_sources; ii++) {
    const Source2D source(unif_rows(rng), unif_cols(rng));
    sources_.push_back(source);
  }

  // Choose a random initial pose.
  pose_ = GridPose2D(unif_rows(rng), unif_cols(rng), unif_angle(rng));
}

// Plan a new trajectory.
bool ExplorerLP::PlanAhead(std::vector<GridPose2D>& trajectory) {
  // Use GUROBI to set up and solve a linear program. Since GUROBI uses
  // exceptions to communicate errors, we enclose the entire planner in
  // a try/catch block.
  try {
    // Create an empty environment and model.
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    // Generate conditional distributions.
    Eigen::MatrixXd pzx;
    Eigen::VectorXd hmz;
    std::vector<unsigned int> trajectory_ids;
    map_.GenerateConditionals(num_samples_, num_steps_, pose_, fov_,
                             pzx, hmz, trajectory_ids);

    // Get vector of objective coefficients.
    const Eigen::VectorXd coefficients = pzx.transpose() * hmz;

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

    // Set up objective function and probability constraint.
    GRBLinExpr objective = 0.0;
    GRBLinExpr constraint = 0.0;
    for (unsigned int ii = 0; ii < variables.size(); ii++) {
      objective += coefficients(ii) * variables[ii];
      constraint += variables[ii];
    }

    model.setObjective(objective, GRB_MAXIMIZE);
    model.addConstr(constraint == 1.0, "probability");

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
}

// Take a step along the given trajectory. Return resulting entropy.
double ExplorerLP::TakeStep(const std::vector<GridPose2D>& trajectory) {
  CHECK(trajectory.size() > 0);

  // Update list of past poses.
  past_poses_.push_back(pose_);

  // Update the current pose.
  pose_ = trajectory[0];

  // Update the map, and return entropy.
  const Sensor2D sensor(pose_, fov_);
  map_.Update(sensor, sources_, true);

  return map_.Entropy();
}

// Compute map entropy.
double ExplorerLP::Entropy() const { return map_.Entropy(); }

// Visualize the current belief state.
void ExplorerLP::Visualize() const {
  glClear(GL_COLOR_BUFFER_BIT);

  // Display each grid cell as a GL_QUAD centered at the appropriate location.
  const Eigen::MatrixXd belief = map_.GetImmutableBelief();

  glBegin(GL_QUADS);
  for (unsigned int ii = 0; ii < map_.GetNumRows(); ii++) {
    for (unsigned int jj = 0; jj < map_.GetNumCols(); jj++) {
      glColor3f(static_cast<GLfloat>(belief(ii, jj)),
                static_cast<GLfloat>(belief(ii, jj)),
                static_cast<GLfloat>(belief(ii, jj)));

      // Bottom left, bottom right, top right, top left.
      glVertex2f(static_cast<GLfloat>(ii), static_cast<GLfloat>(jj));
      glVertex2f(static_cast<GLfloat>(ii) + 1.0, static_cast<GLfloat>(jj));
      glVertex2f(static_cast<GLfloat>(ii) + 1.0, static_cast<GLfloat>(jj) + 1.0);
      glVertex2f(static_cast<GLfloat>(ii), static_cast<GLfloat>(jj) + 1.0);
    }
  }
  glEnd();

  const GLfloat robot_x = static_cast<GLfloat>(pose_.GetX());
  const GLfloat robot_y = static_cast<GLfloat>(pose_.GetY());
  const GLfloat robot_a = static_cast<GLfloat>(pose_.GetAngle());
  const unsigned int kNumVertices = 100;

  // Display the field of view as a triangle fan.
  const GLfloat kFovRadius =
    sqrt(static_cast<GLfloat>(map_.GetNumRows() * map_.GetNumRows() +
                            map_.GetNumCols() * map_.GetNumCols()));

  glBegin(GL_TRIANGLE_FAN);
  glColor4f(0.0, 0.2, 0.8, 0.25);
  glVertex2f(robot_x, robot_y);
  for (unsigned int ii = 0; ii <= kNumVertices; ii++) {
    const GLfloat angle = robot_a + fov_ *
      (-0.5 + static_cast<GLfloat>(ii) / static_cast<GLfloat>(kNumVertices));
    glVertex2f(robot_x + kFovRadius * cos(angle),
               robot_y + kFovRadius * sin(angle));
  }
  glEnd();

  // Display a circle at the robot's current position. No circle primitive, so
  // use a polygon with a bunch of vertices.
  const GLfloat kRobotRadius = 0.5;

  glBegin(GL_POLYGON);
  glColor4f(0.0, 0.8, 0.2, 0.75);
  for (unsigned int ii = 0; ii < kNumVertices; ii++) {
    const GLfloat angle =
      2.0 * M_PI * static_cast<GLfloat>(ii) / static_cast<GLfloat>(kNumVertices);
    glVertex2f(robot_x + kRobotRadius * cos(angle),
               robot_y + kRobotRadius * sin(angle));
  }
  glEnd();
}

} // namespace radiation
