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
// Defines a problem setup, parameterized by a number of sources, number of
// steps in each trajectory, sensor field of view, and number of samples to
// estimate distributions from.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIATION_PROBLEM_2D_H
#define RADIATION_PROBLEM_2D_H

#include "source_2d.h"
#include "sensor_2d.h"
#include "grid_map_2d.h"
#include "grid_pose_2d.h"
#include "movement_2d.h"
#include "encoding.h"

#include <Eigen/Core>

namespace radiation {

class Problem2D {
 public:
  Problem2D(unsigned int num_sources, unsigned int num_steps, double fov,
            unsigned int num_samples, );
  ~Problem2D();


 private:
  // Solve least squares problem to update belief state.
  void SolveLeastSquares();

  // Belief state.
  Eigen::VectorXd belief_;

  // Problem parameters.
  const unsigned int num_rows_;
  const unsigned int num_cols_;
  const unsigned int num_sources_;

  // Regularizer for belief update. Enforeces consistency across all voxels.
  const double regularizer_;

  // List of viewed indices, where each index is represented as a tuple, and
  // corresponding measurements.
  std::vector<std::vector<std::tuple<unsigned int, unsigned int>>> viewed_;
  std::vector<unsigned int> measurements_;
}; // class GridMap2D

} // namespace radiation

#endif
