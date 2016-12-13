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

#include <explorer_socp.h>
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

DEFINE_string(rw_output_file, "rw_out.csv",
              "Name of file to save random walk results in.");
DEFINE_string(lp_output_file, "lp_out.csv",
              "Name of file to save LP results in.");
DEFINE_string(socp_output_file, "socp_out.csv",
              "Name of file to save SOCP results in.");
DEFINE_bool(verbose, false, "Output iteration number and final entropies.");
DEFINE_int32(num_trials, 1000, "Number of random trials to run.");
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
DEFINE_double(epsilon, 0.1, "Noise ball radius for SOCP explorer.");

using namespace radiation;

// Set everything up and go!
int main(int argc, char** argv) {
  // Set up logging.
  google::InitGoogleLogging(argv[0]);

  // Parse flags.
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Make sure we can open the output files.
  std::ofstream rw_file(FLAGS_rw_output_file);
  if (!rw_file.is_open()) {
    VLOG(1) << "Unable to open file " + FLAGS_rw_output_file;
    return 1;
  }

  std::ofstream lp_file(FLAGS_lp_output_file);
  if (!lp_file.is_open()) {
    VLOG(1) << "Unable to open file " + FLAGS_lp_output_file;
    return 1;
  }

  std::ofstream socp_file(FLAGS_socp_output_file);
  if (!socp_file.is_open()) {
    VLOG(1) << "Unable to open file " + FLAGS_socp_output_file;
    return 1;
  }

  // Set static variables.
  GridPose2D::SetNumRows(FLAGS_num_rows);
  GridPose2D::SetNumCols(FLAGS_num_cols);
  Movement2D::SetAngularStep(FLAGS_angular_step);

  // Set up a random number generator.
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<unsigned int> unif_rows(0, FLAGS_num_rows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, FLAGS_num_cols - 1);
  std::uniform_real_distribution<double> unif_angle(0.0, 2.0 * M_PI);

  // Only plan ahead if we haven't exceeded the step count.
  for (unsigned int ii = 0; ii < FLAGS_num_trials; ii++) {
    // Create a random set of sources.
    std::vector<Source2D> sources;
    for (unsigned int jj = 0; jj < FLAGS_num_sources; jj++) {
      const Source2D source(unif_rows(rng), unif_cols(rng));
      sources.push_back(source);
    }

    // Create a random initial pose.
    const GridPose2D initial_pose(unif_rows(rng), unif_cols(rng),
                                  unif_angle(rng));

    // Set up explorers.
    ExplorerRW rw(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                  FLAGS_regularizer, FLAGS_fov, sources, initial_pose);
    ExplorerLP lp(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                  FLAGS_regularizer, FLAGS_num_steps, FLAGS_fov,
                  FLAGS_num_samples, sources, initial_pose);
    ExplorerSOCP socp(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                      FLAGS_regularizer, FLAGS_num_steps, FLAGS_fov,
                      FLAGS_num_samples, FLAGS_epsilon, sources, initial_pose);

    // Write initial entropies to file (these will be the same).
    double rw_entropy = rw.Entropy();
    rw_file << rw_entropy;

    double lp_entropy = lp.Entropy();
    lp_file << lp_entropy;

    double socp_entropy = socp.Entropy();
    socp_file << socp_entropy;

    // Run for the specified number of steps.
    std::vector<GridPose2D> lp_trajectory, rw_trajectory, socp_trajectory;
    for (unsigned int jj = 0; jj < FLAGS_num_iterations; jj++) {
      // Plan ahead.
      CHECK(lp.PlanAhead(lp_trajectory));
      CHECK(rw.PlanAhead(rw_trajectory));
      CHECK(socp.PlanAhead(socp_trajectory));

      // Take a step.
      lp_entropy = lp.TakeStep(lp_trajectory);
      rw_entropy = rw.TakeStep(rw_trajectory);
      socp_entropy = socp.TakeStep(socp_trajectory);

      // Write to files.
      rw_file << ", " << rw_entropy;
      lp_file << ", " << lp_entropy;
      socp_file << ", " << socp_entropy;
    }

    // Write newlines.
    rw_file << "\n";
    lp_file << "\n";
    socp_file << "\n";

    // Maybe output to screen.
    if (FLAGS_verbose && ii % 100 == 0) {
      std::cout << "Trial " << ii << ": RW = " << rw_entropy <<
        ", LP = " << lp_entropy << ", SOCP = " << socp_entropy <<
        std::endl << std::flush;
    }
  }

  // Close files.
  rw_file.close();
  lp_file.close();
  socp_file.close();
  return 0;
}
