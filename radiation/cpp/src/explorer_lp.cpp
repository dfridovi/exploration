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

#include <glog/logging.h>
#include <gurobi_c++.h>
#include <random>
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
    fov_(fov) {
  // Set up a random number generator.
  std::random_device rd;
  std::default_random_generator rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, num_rows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, num_cols - 1);
  std::uniform_real_distibution<double> unif_angle(0.0, 2.0 * M_PI);

  // Choose random sources.
  for (unsigned int ii = 0; ii < num_sources; ii++) {
    const Source2D source(unif_rows(rng), unif_cols(rng));
    sources_.push_back(source);
  }

  // Choose a random initial pose.
  pose_ = GridPose2D(unif_rows(rng), unif_cols(rng), unif_angle(rng));
}

// Plan a new trajectory.
bool ExplorerLP::PlanAhead(std::vector<GridPose2D>& trajectory) const {
  // TODO!
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

// Visualize the current belief state.
void ExplorerLP::Visualize(const std::string& title) const {
  // TODO!
}

} // namespace radiation
