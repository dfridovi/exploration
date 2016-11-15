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

#include "grid_map_2d.h"

#include <glog/logging.h>
#include <algorithm>
#include <math.h>

namespace radiation {

  GridMap2D::~GridMap2D() {}
  GridMap2D::GridMap2D(unsigned int num_rows, unsigned int num_cols,
                       unsigned int num_sources, double regularizer)
    : num_rows_(num_rows), num_cols_(num_cols_),
      num_sources_(num_sources), regularizer_(regularizer),
      rng_(rd_) {

    // Initialize belief matrix to be uniform.
    belief_ = Eigen::MatrixXd::Ones(num_rows_, num_cols_) * num_sources_;
    belief_ /= num_rows_ * num_cols_;
  }

  // Generate random sources according to the current belief state.
  bool GridMap2D::GenerateSources(std::vector<Source2D>& sources) const {
    const double total_belief = belief_.sum();
    std::uniform_real_distribution unif(0.0, 1.0);

    // Choose 'num_sources_' random numbers in [0, 1], which will be sorted
    // and treated as evaluations of the CDF. Since they are uniform, the points
    // at which they occur ar distributed according to the current belief.
    std::vector<double> cdf_evals;
    for (size_t ii = 0; ii < num_sources_; ii++)
      cdf_evals.push_back(unif(rng_));

    std::sort(cdf_evals.begin(), cdf_evals.end());

    // Walk the current belief distribution until we get to each 'cdf_eval'
    // and record which voxel we are in.
    sources.clear();
    unsigned int current_index = 0;
    double current_cdf = 0.0;

    for (size_t ii = 0; ii < num_rows_; ii++) {
      for (size_t jj = 0; jj < num_cols_; jj++) {
        current_cdf += belief_(ii, jj) / total_belief;

        // Check if we just passed the next 'cdf_eval'.
        while (current_cdf >= cdf_evals[current_index]) {
          // Generate a new source here.
          sources.push_back(Source2D(ii, jj));

          // Return if we have enough sources.
          if (sources.size() == num_sources_)
            return true;

          // Increment 'current_index'.
          current_index++;
        }
      }
    }

    // Should never get here.
    return false;
  }

  // Generate conditional distribution [P_{Z|X}] and entropy vector [h_{M|Z}],
  // where the (i, j)-entry of [P_{Z|X}] is the normalized frequency of
  // observing measurement i given trajectory j, and the i-entry of [h_{M|Z}]
  // is the entropy of M given measurement i, starting from the given pose.
  void GridMap2D::GenerateConditionals(
     unsigned int num_samples, const GridPose2D& pose,
     Eigen::MatrixXd& pzx, Eigen::VectorXd& hmz,
     std::vector<unsigned int>& trajectory_ids) const {
    // TODO!
  }

  // Take a measurement from the given sensor and update belief accordingly.
  bool GridMap2D::Update(const Sensor2D& sensor,
                         const std::vector<Source2D>& sources,
                         bool solve = true) {
    const unsigned int measurement = sensor.Sense(sources);
    CHECK(measurement <= num_sources_);

    // Identify all voxels in range and store.
    std::vector< std::tuple<unsigned int, unsigned int> > voxels;
    for (size_t ii = 0; ii < num_rows_; ii++) {
      for (size_t jj = 0; jj < num_cols_; jj++) {
        if (sensor.VoxelInView(ii, jj))
          voxels.push_back(std::make_tuple(ii, jj));
      }
    }

    viewed_.push_back(voxels);
    measurements_.push_back(measurement);

    // Maybe solve.
    if (solve)
      return SolveLeastSquares();

    return true;
  }

  // Compute entropy.
  double GridMap2D::Entropy() const {
    double entropy = 0.0;

    for (size_t ii = 0; ii < num_rows_; ii++) {
      for (size_t jj = 0; jj < num_cols_; jj++) {
        const double p = belief_(ii, jj);
        const double bernoulli_entropy = (p > 1e-8 && p < 1.0 - 1e-8) ?
          (-p * log(p) - (1.0 - p) * log(1.0 - p)) : 0.0;

        entropy += bernoulli_entropy;
      }
    }

    return entropy;
  }

  // Solve least squares problem to update belief state.
  bool GridMap2D::SolveLeastSquares() {
    // TODO!
  }

} // namespace radiation
