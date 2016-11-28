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

#include <explorer_lp.h>
#include <movement_2d.h>
#include <grid_pose_2d.h>

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <iostream>
#include <math.h>

DEFINE_int32(num_iterations, 10, "Number of iterations to run exploration.");
DEFINE_int32(num_rows, 5, "Number of rows in the grid.");
DEFINE_int32(num_cols, 5, "Number of columns in the grid.");
DEFINE_int32(num_sources, 2, "Number of sources on the grid.");
DEFINE_int32(num_steps, 3, "Number of steps in each trajectory.");
DEFINE_int32(num_samples, 10000,
              "Number of samples used to approximate distributions.");
DEFINE_double(angular_step, 0.33 * M_PI, "Angular step size.");
DEFINE_double(fov, 0.5 * M_PI, "Sensor field of view.");
DEFINE_double(regularizer, 1.0, "Regularization parameter for belief update.");

using namespace radiation;

int main(int argc, char** argv) {
  // Set up logging.
  google::InitGoogleLogging(argv[0]);

  // Parse flags.
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Set static variables.
  GridPose2D::SetNumRows(FLAGS_num_rows);
  GridPose2D::SetNumCols(FLAGS_num_cols);
  Movement2D::SetAngularStep(FLAGS_angular_step);

  // Create an explorer.
  ExplorerLP explorer(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                      FLAGS_regularizer, FLAGS_num_steps, FLAGS_fov,
                      FLAGS_num_samples);

  // Explore for the specified number of iterations.
  for (unsigned int ii = 0; ii < FLAGS_num_iterations; ii++) {
    std::vector<GridPose2D> trajectory;
    if (!explorer.PlanAhead(trajectory)) {
      VLOG(1) << "Explorer encountered an error. Skipping this iteration.";
      continue;
    }

    // Take a step.
    const double entropy = explorer.TakeStep(trajectory);
    std::printf("Entropy after step %u is %f.\n", ii, entropy);

    // Visualize.
    explorer.Visualize();
  }

  return 0;
}
