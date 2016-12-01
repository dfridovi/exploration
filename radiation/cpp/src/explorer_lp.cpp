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
  : Explorer2D(num_rows, num_cols, num_sources, regularizer, fov),
    num_steps_(num_steps),
    num_samples_(num_samples) {}

// Plan a new trajectory.
bool ExplorerLP::PlanAhead(std::vector<GridPose2D>& trajectory) {
  // Generate conditional entropy vector.
  Eigen::VectorXd hzx;
  std::vector<unsigned int> trajectory_ids;
  map_.GenerateEntropyVector(num_samples_, num_steps_, pose_, fov_,
                             hzx, trajectory_ids);
  CHECK(hzx.rows() == trajectory_ids.size());

  // Compute the arg max of this conditional entropy vector.
  double max_value = -1.0;
  unsigned int trajectory_id = 0;
  for (unsigned int ii = 0; ii < hzx.rows(); ii++) {
    if (hzx(ii) > max_value) {
      max_value = hzx(ii);
      trajectory_id = trajectory_ids[ii];
    }
  }

  // Check that we found a valid trajectory (with non-negative entropy).
  if (max_value < 0.0) {
    VLOG(1) << "Could not find a positive conditional entropy trajectory.";
    return false;
  }

  // Decode this trajectory id.
  trajectory.clear();
  DecodeTrajectory(trajectory_id, num_steps_, pose_, trajectory);
  return true;
}

} // namespace radiation
