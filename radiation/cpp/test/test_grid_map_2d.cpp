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
#include <math.h>

namespace radiation {

// Test ethat we can detect sources randomly located across the grid.
TEST(GridMap2D, TestEmptyMap) {
  const unsigned int kNumRows = 10;
  const unsigned int kNumCols = 10;
  const unsigned int kNumSources = 1;
  const double kRegularizer = 0.0;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);

  // Make random number generators.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);
  std::uniform_real_distribution<double> unif_angle(0.0, 0.25 * M_PI);

  // Always place the sensor at the same pose.
  const double kAngle = 0.25 * M_PI;
  const GridPose2D sensor_pose(0.0, 0.0, kAngle);

  // Generate a buch of random sources.
  std::vector<Source2D> sources;
  for (size_t jj = 0; jj < kNumSources; jj++)
    sources.push_back(Source2D(unif_rows(rng), unif_cols(rng)));

  // Generate a random field of view..
  const double fov = unif_angle(rng);
  const Sensor2D sensor(sensor_pose, fov);

  // Count the sources in view, using the fact that the FOV is a cone
  // looking into the first quadrant.
  unsigned int total_count = 0;
  for (const auto& source : sources) {
    const double angle = atan2(source.GetY(), source.GetX());
    if ((angle >= kAngle - 0.5 * fov) && (angle <= kAngle + 0.5 * fov)) {
      EXPECT_TRUE(sensor.SourceInView(source));
      total_count++;
    } else {
      EXPECT_FALSE(sensor.SourceInView(source));
    }
  }

  // Check that the sensor measurement matches 'total_count'.
  EXPECT_EQ(sensor.Sense(sources), total_count);
}

} // namespace radiation
