/*
 * Copyright (c) 2016, The Regents of the University of California (Regents).
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
// Unit tests for the GridMap2D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <grid_map_2d.h>
#include <sensor_2d.h>
#include <source_2d.h>

#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <iostream>
#include <limits>
#include <math.h>

namespace radiation {

// Test that we can detect sources randomly located across the grid.
TEST(GridMap2D, TestEmptyMap) {
  const unsigned int kNumRows = 10;
  const unsigned int kNumCols = 10;
  const unsigned int kNumSources = 1;
  const double kRegularizer = 0.0;
  const unsigned int kNumUpdates = 10;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);

  // Make random number generators.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);

  // Place an omnidirectional sensor at the origin.
  const double kAngle = 0.0;
  const double kFov = 2.0 * M_PI;
  const GridPose2D sensor_pose(0.0, 0.0, kAngle);
  const Sensor2D sensor(sensor_pose, kFov);

  // Update a bunch of times, and check that the maps' 'belief' has
  // converged to zero.
  std::vector<Source2D> empty_sources;
  GridMap2D map(kNumRows, kNumCols, kNumSources, kRegularizer);

  double belief_norm = map.GetImmutableBelief().norm();
  for (unsigned int ii = 0; ii < kNumUpdates; ii++) {
    EXPECT_TRUE(map.Update(sensor, empty_sources, true));

    // Make sure norm of 'belief' has decreased.
    const double updated_norm = map.GetImmutableBelief().norm();
    EXPECT_LE(updated_norm, belief_norm);
    belief_norm = updated_norm;
  }

  const Eigen::MatrixXd belief = map.GetImmutableBelief();
  for (unsigned int ii = 0; ii < kNumRows; ii++)
    for (unsigned int jj = 0; jj < kNumCols; jj++)
      EXPECT_NEAR(belief(ii, jj), 0.0, 1e-4);
}

// Test that we can detect sources randomly located across the grid.
TEST(GridMap2D, TestConvergenceSingleSource) {
  const unsigned int kNumRows = 5;
  const unsigned int kNumCols = 5;
  const unsigned int kNumSources = 1;
  const double kRegularizer = 1.0;
  const unsigned int kNumUpdates = 100;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);

  // Make random number generators.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);
  std::uniform_real_distribution<double> unif_angle(0.0, 2.0 * M_PI);

  // Choose random source locations.
  std::vector<Source2D> sources;
  for (unsigned int ii = 0; ii < kNumSources; ii++)
    sources.push_back(Source2D(unif_rows(rng), unif_cols(rng)));

  // Create a new map.
  GridMap2D map(kNumRows, kNumCols, kNumSources, kRegularizer);

  // Iterate the specified number of times. Each time, choose a random sensor
  // pose and take a measurement. Update the map and repeat.
  const double kFov = 0.2 * M_PI;
  double entropy = map.Entropy();
  for (unsigned int ii = 0; ii < kNumUpdates; ii++) {
    const GridPose2D pose(unif_rows(rng), unif_cols(rng), unif_angle(rng));
    const Sensor2D sensor(pose, kFov);

    EXPECT_TRUE(map.Update(sensor, sources, true));

    // Make sure entropy has not increased by much.
    const double updated_entropy = map.Entropy();
    EXPECT_LE(updated_entropy, 2.0 * entropy);
    entropy = updated_entropy;
  }

  // Check that belief has converged to the truth.
  const Eigen::MatrixXd belief = map.GetImmutableBelief();
  for (unsigned int ii = 0; ii < kNumRows; ii++) {
    for (unsigned int jj = 0; jj < kNumCols; jj++) {
      // If there is a source here, then make sure 'belief' has found it.
      bool has_source = false;
      for (const auto& source : sources) {
        if (source.GetIndexX() == ii && source.GetIndexY() == jj) {
          has_source = true;
          break;
        }
      }

      if (has_source)
        EXPECT_GE(belief(ii, jj), 1.0 - 1e-4);
      else
        EXPECT_LE(belief(ii, jj), 1e-4);
    }
  }
}

// Test that we can detect sources randomly located across the grid.
TEST(GridMap2D, TestConvergenceMultipleSources) {
  const unsigned int kNumRows = 5;
  const unsigned int kNumCols = 5;
  const unsigned int kNumSources = 2;
  const double kRegularizer = 1.0;
  const unsigned int kNumUpdates = 200;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);

  // Make random number generators.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);
  std::uniform_real_distribution<double> unif_angle(0.0, 2.0 * M_PI);

  // Choose random source locations.
  std::vector<Source2D> sources;
  for (unsigned int ii = 0; ii < kNumSources; ii++)
    sources.push_back(Source2D(unif_rows(rng), unif_cols(rng)));

  // Create a new map.
  GridMap2D map(kNumRows, kNumCols, kNumSources, kRegularizer);

  // Iterate the specified number of times. Each time, choose a random sensor
  // pose and take a measurement. Update the map and repeat.
  const double kFov = 0.2 * M_PI;
  double entropy = map.Entropy();
  for (unsigned int ii = 0; ii < kNumUpdates; ii++) {
    const GridPose2D pose(unif_rows(rng), unif_cols(rng), unif_angle(rng));
    const Sensor2D sensor(pose, kFov);

    EXPECT_TRUE(map.Update(sensor, sources, true));

    // Make sure entropy has not increased by much.
    const double updated_entropy = map.Entropy();
    EXPECT_LE(updated_entropy, 2.0 * entropy);
    entropy = updated_entropy;
  }

  // Check that belief has converged to the truth.
  const Eigen::MatrixXd belief = map.GetImmutableBelief();
  for (unsigned int ii = 0; ii < kNumRows; ii++) {
    for (unsigned int jj = 0; jj < kNumCols; jj++) {
      // If there is a source here, then make sure 'belief' has found it.
      bool has_source = false;
      for (const auto& source : sources) {
        if (source.GetIndexX() == ii && source.GetIndexY() == jj) {
          has_source = true;
          break;
        }
      }

      if (has_source)
        EXPECT_GE(belief(ii, jj), 1.0 - 1e-4);
      else
        EXPECT_LE(belief(ii, jj), 1e-4);
    }
  }
}

// Test that we can detect sources randomly located across the grid.
TEST(GridMap2D, TestDistributionConvergence) {
  const unsigned int kNumRows = 5;
  const unsigned int kNumCols = 5;
  const unsigned int kNumSources = 1;
  const unsigned int kNumSteps = 1;
  const unsigned int kNumSamples = 100000;
  const double kAngularStep = 0.25 * M_PI;
  const double kFov = 0.2 * M_PI;
  const double kPrecision = 0.1;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);
  Movement2D::SetAngularStep(kAngularStep);

  // Make random number generators.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);
  std::uniform_int_distribution<unsigned int> unif_angle(0.0, 2.0 * M_PI);

  // Create a new map.
  GridMap2D map(kNumRows, kNumCols, kNumSources, 0.0 /* regularizer */);

  // Create a random initial pose.
  const GridPose2D pose(unif_rows(rng), unif_cols(rng), unif_angle(rng));

  // Generate conditionals twice.
  Eigen::MatrixXd pzx1, pzx2;
  Eigen::VectorXd hmz1, hmz2;
  std::vector<unsigned int> trajectory_ids1, trajectory_ids2;

  map.GenerateConditionals(kNumSamples, kNumSteps, pose, kFov,
                           pzx1, hmz1, trajectory_ids1);
  map.GenerateConditionals(kNumSamples, kNumSteps, pose, kFov,
                           pzx2, hmz2, trajectory_ids2);

  // Check that trajectory ids match. This is a simple check since
  // they are generated from an _ordered_ map.
  EXPECT_EQ(trajectory_ids1.size(), trajectory_ids2.size());

  if (trajectory_ids1.size() == trajectory_ids2.size()) {
    ASSERT_EQ(pzx1.rows(), pzx2.rows());
    ASSERT_EQ(pzx1.cols(), pzx2.cols());
    ASSERT_EQ(hmz1.rows(), hmz2.rows());
    ASSERT_EQ(hmz1.cols(), hmz2.cols());
    ASSERT_EQ(pzx1.rows(), hmz1.rows());

    bool ids_match = true;
    for (unsigned int ii = 0; ii < trajectory_ids1.size(); ii++) {
      EXPECT_EQ(trajectory_ids1[ii], trajectory_ids2[ii]);
      ids_match &= (trajectory_ids1[ii] == trajectory_ids2[ii]);
    }

    if (ids_match) {
      EXPECT_TRUE(pzx1.isApprox(pzx2, kPrecision));
      EXPECT_TRUE(hmz1.isApprox(hmz2, kPrecision));

      const Eigen::VectorXd c1 = pzx1.transpose() * hmz1;
      const Eigen::VectorXd c2 = pzx2.transpose() * hmz2;
      EXPECT_TRUE(c1.isApprox(c2, kPrecision));
    }
  }
}

} // namespace radiation
