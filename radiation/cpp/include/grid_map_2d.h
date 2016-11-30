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
// Defines a 2D grid map.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIATION_GRID_MAP_2D_H
#define RADIATION_GRID_MAP_2D_H

#include <source_2d.h>
#include <sensor_2d.h>
#include <grid_pose_2d.h>

#include <Eigen/Core>

#include <random>
#include <vector>

namespace radiation {

class GridMap2D {
 public:
  GridMap2D(unsigned int num_rows, unsigned int num_cols,
            unsigned int num_sources, double regularizer);
  ~GridMap2D();

  // Getters.
  unsigned int GetNumRows() const;
  unsigned int GetNumCols() const;

  // Generate random sources according to the current belief state.
  bool GenerateSources(std::vector<Source2D>& sources);

  // Generate entropy vector [h_{Z|X}], where the i-entry of [h_{Z|X}]
  // is the entropy of Z given trajectory X = i, starting from the given pose.
  void GenerateEntropyVector(unsigned int num_samples, unsigned int num_steps,
                            const GridPose2D& pose, double sensor_fov,
                            Eigen::VectorXd& hzx,
                            std::vector<unsigned int>& trajectory_ids);

  // Take a measurement from the given sensor and update belief accordingly.
  bool Update(const Sensor2D& sensor,
              const std::vector<Source2D>& sources, bool solve = true);

  // Compute entropy.
  double Entropy() const;

  // Get a reference to immutable 'belief'.
  const Eigen::MatrixXd& GetImmutableBelief() const;

 private:
  // Solve least squares problem to update belief state.
  bool SolveLeastSquares();

  // Belief state.
  Eigen::MatrixXd belief_;

  // Problem parameters.
  const unsigned int num_rows_;
  const unsigned int num_cols_;
  const unsigned int num_sources_;

  // Regularizer for belief update. Enforeces consistency across all voxels.
  const double regularizer_;

  // List of viewed indices, where each index is given in column-major form, and
  // corresponding measurements.
  std::vector< std::vector<unsigned int> > viewed_;
  std::vector<unsigned int> measurements_;

  // Random number generator.
  std::random_device rd_;
  std::default_random_engine rng_;
}; // class GridMap2D

} // namespace radiation

#endif
