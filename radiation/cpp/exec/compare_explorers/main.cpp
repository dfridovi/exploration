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
#include <explorer_rw.h>
#include <movement_2d.h>
#include <grid_pose_2d.h>

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <iostream>
#include <fstream>
#include <math.h>

#ifdef SYSTEM_OSX
#include <GLUT/glut.h>
#endif

#ifdef SYSTEM_LINUX
#include <GL/glut.h>
#endif

DEFINE_string(output_file, "out.csv",
              "Name of file to save results in. Column 1 is RW, 2 is LP.");
DEFINE_bool(verbose, false, "Output iteration number and final entropies.");
DEFINE_int32(num_trials, 10000, "Number of random trials to run.");
DEFINE_int32(num_iterations, 10, "Number of iterations to run exploration.");
DEFINE_int32(num_rows, 5, "Number of rows in the grid.");
DEFINE_int32(num_cols, 5, "Number of columns in the grid.");
DEFINE_int32(num_sources, 2, "Number of sources on the grid.");
DEFINE_int32(num_steps, 2, "Number of steps in each trajectory.");
DEFINE_int32(num_samples, 20000,
              "Number of samples used to approximate distributions.");
DEFINE_double(angular_step, 0.07 * M_PI, "Angular step size.");
DEFINE_double(fov, 0.1 * M_PI, "Sensor field of view.");
DEFINE_double(regularizer, 1.0, "Regularization parameter for belief update.");

using namespace radiation;

// Set everything up and go!
int main(int argc, char** argv) {
  // Set up logging.
  google::InitGoogleLogging(argv[0]);

  // Parse flags.
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Make sure we can open the output file.
  std::ofstream file(FLAGS_output_file);
  if (!file.is_open()) {
    VLOG(1) << "Unable to open file " + FLAGS_output_file;
    return 1;
  }

  // Set static variables.
  GridPose2D::SetNumRows(FLAGS_num_rows);
  GridPose2D::SetNumCols(FLAGS_num_cols);
  Movement2D::SetAngularStep(FLAGS_angular_step);

  // Only plan ahead if we haven't exceeded the step count.
  for (unsigned int ii = 0; ii < FLAGS_num_trials; ii++) {
    // Set up explorers.
    ExplorerLP lp(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                  FLAGS_regularizer, FLAGS_num_steps, FLAGS_fov,
                  FLAGS_num_samples);
    ExplorerRW rw(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                  FLAGS_regularizer, FLAGS_fov);

    double rw_entropy = rw.Entropy();
    double lp_entropy = lp.Entropy();

    // Run for the specified number of steps.
    for (unsigned int jj = 0; jj < FLAGS_num_iterations; jj++) {
      // Plan ahead.
      std::vector<GridPose2D> lp_trajectory, rw_trajectory;
      CHECK(lp.PlanAhead(lp_trajectory));
      CHECK(rw.PlanAhead(rw_trajectory));

      // Take a step.
      lp_entropy = lp.TakeStep(lp_trajectory);
      rw_entropy = rw.TakeStep(rw_trajectory);
    }

    // Write to file.
    file << rw_entropy << ", " << lp_entropy << "\n";

    // Maybe output to screen.
    if (FLAGS_verbose) {
      std::cout << "Trial " << ii << ": RW = " << rw_entropy <<
        ", LP = " << lp_entropy << std::endl << std::flush;
    }
  }

  // Close file.
  file.close();
  return 0;
}
