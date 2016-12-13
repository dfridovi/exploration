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
// Unit tests for the ExplorerSOCP class.
//
///////////////////////////////////////////////////////////////////////////////

#include <grid_pose_2d.h>
#include <explorer_socp.h>
#include <explorer_rw.h>

#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <iostream>
#include <limits>
#include <math.h>

namespace radiation {

// Test that the SOCP explorer does better than a random walk. Since each
// explorer initializes with a random sample configuration, we test
// performance after a specified number of iterations, averaged over a
// large number of trials. Performance is measured as total entropy
// over all iterations, averaged over trials.
TEST(ExplorerSOCP, TestVsRandomWalk) {
  const unsigned int kNumRows = 1;
  const unsigned int kNumCols = 20;
  const unsigned int kNumSources = 1;
  const unsigned int kNumSteps = 3;
  const unsigned int kNumSamples = 10000;
  const double kRegularizer = 1.0;
  const unsigned int kNumIterations = 3;
  const unsigned int kNumTrials = 10;
  const double kAngularStep = 0.07 * M_PI;
  const double kFov = 0.1 * M_PI;
  const double kEpsilon = 0.1;

  // Set static variables.
  GridPose2D::SetNumRows(kNumRows);
  GridPose2D::SetNumCols(kNumCols);
  Movement2D::SetAngularStep(kAngularStep);

  // Create explorers.
  ExplorerSOCP socp(kNumRows, kNumCols, kNumSources, kRegularizer,
                    kNumSteps, kFov, kNumSamples, kEpsilon);
  ExplorerRW rw(kNumRows, kNumCols, kNumSources, kRegularizer, kFov);

  // Average over a bunch of trials.
  double socp_entropy = 0.0;
  double rw_entropy = 0.0;
  std::vector<GridPose2D> trajectory;
  for (unsigned int ii = 0; ii < kNumTrials; ii++) {
    for (unsigned int jj = 0; jj < kNumIterations; jj++) {
      // SOCP explorer.
      EXPECT_TRUE(socp.PlanAhead(trajectory));
      socp_entropy += socp.TakeStep(trajectory);

      // RW explorer.
      EXPECT_TRUE(rw.PlanAhead(trajectory));
      rw_entropy += rw.TakeStep(trajectory);
    }
  }

  EXPECT_LE(socp_entropy / static_cast<double>(kNumTrials),
            rw_entropy / static_cast<double>(kNumTrials));
}

} // namespace radiation
