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
 * Authors: David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 */

///////////////////////////////////////////////////////////////////////////////
//
// This file defines cost functors that can be used in conjunction with Google
// Ceres solver to solve non-linear least-squares problems. Each functor must
// define a public templated operator() method that takes in arguments const T*
// const INPUT_VARIABLE, and T* OUTPUT_RESIDUAL, and returns a bool.
// 'INPUT_VARIABLE' will be the optimization variable, and 'OUTPUT_RESIDUAL'
// will be the cost associated with the input.
//
// One can define more specific cost functions by adding other structure to the
// functor, e.g. by passing in other parameters of the cost function to the
// functor's constructor.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIATION_COST_FUNCTORS_H
#define RADIATION_COST_FUNCTORS_H

#include <ceres/ceres.h>
#include <Eigen/Core>
#include <glog/logging.h>
#include <vector>
#include <math.h>

namespace radiation {

// Belief error measures the difference between the expected sensor measurement
// and the true measurement, assuming the given belief vector of probabilities
// that each voxel is occupied.
struct BeliefError {
  // Inputs: pointer to list of voxel indices that are in view, and measurement.
  // Optimization variables are the probability that each voxel contains a
  // source.
  const std::vector<unsigned int>* voxels_;
  const unsigned int measurement_;

  BeliefError(const std::vector<unsigned int>* voxels,
              unsigned int measurement)
    : voxels_(voxels), measurement_(measurement) {
    CHECK_NOTNULL(voxels);
  }

  template <typename T>
  bool operator()(const T* const belief, T* expected_error) const {
    // Compute the difference between the expected measurement and the
    // actual measurement.
    *expected_error = -static_cast<T>(measurement_);

    for (const auto& voxel_id : *voxels_) {
      *expected_error += belief[voxel_id];
    }

    return true;
  }

  // Factory method.
  static ceres::CostFunction* Create(const std::vector<unsigned int>* voxels,
                                     const unsigned int measurement) {
    // Only a single residual.
    static const int kNumResiduals = 1;

    // Number of parameters is the number of grid cells.
    // MAKE SURE TO CHANGE THIS IF THE NUMBER OF GRID CELLS CHANGES.
    static const int kNumParameters = 25;

    return new ceres::AutoDiffCostFunction<BeliefError,
      kNumResiduals,
      kNumParameters>(new BeliefError(voxels, measurement));
  }
};  // struct BeliefError

// Belief regularization enforces consistency across the entire belief vector
// by measuring the difference between the expected number of sources and the
// total number that should be present.
struct BeliefRegularization {
  // Inputs: true number of sources and regularization tradeoff parameter.
  // Optimization variables are the probability that each voxel contains a
  // source.
  const unsigned int num_sources_;
  const double regularizer_;

  BeliefRegularization(unsigned int num_sources, double regularizer)
  : num_sources_(num_sources), regularizer_(sqrt(regularizer)) {}

  template <typename T>
  bool operator()(const T* const belief, T* expected_error) const {
    // Compute the difference between the expected measurement and the
    // actual measurement.
    *expected_error = -static_cast<T>(num_sources_);

    // Set this equal to 25 for now.
    const unsigned int kNumParameters = 25;

    for (size_t ii = 0; ii < kNumParameters; ii++) {
      *expected_error += belief[ii];
    }

    *expected_error *= static_cast<T>(regularizer_);

    return true;
  }

  // Factory method.
  static ceres::CostFunction* Create(unsigned int num_sources,
                                     double regularizer) {
    // Only a single residual.
    static const int kNumResiduals = 1;

    // Number of parameters is the number of grid cells.
    // MAKE SURE TO CHANGE THIS IF THE NUMBER OF GRID CELLS CHANGES.
    static const int kNumParameters = 25;

    return new ceres::AutoDiffCostFunction<BeliefRegularization,
      kNumResiduals,
      kNumParameters>(new BeliefRegularization(num_sources, regularizer));
  }
};  // struct BeliefRegularization

}  // namespace radiation

#endif
