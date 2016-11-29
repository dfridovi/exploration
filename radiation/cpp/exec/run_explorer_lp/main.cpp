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

#include <GLUT/glut.h>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <iostream>
#include <math.h>

DEFINE_int32(refresh_rate, 1000, "Refresh rate in milliseconds.");
DEFINE_bool(iterate_forever, false, "Iterate ad inifinitum?");
DEFINE_int32(num_iterations, 10, "Number of iterations to run exploration.");
DEFINE_int32(num_rows, 5, "Number of rows in the grid.");
DEFINE_int32(num_cols, 5, "Number of columns in the grid.");
DEFINE_int32(num_sources, 2, "Number of sources on the grid.");
DEFINE_int32(num_steps, 4, "Number of steps in each trajectory.");
DEFINE_int32(num_samples, 20000,
              "Number of samples used to approximate distributions.");
DEFINE_double(angular_step, 0.07 * M_PI, "Angular step size.");
DEFINE_double(fov, 0.1 * M_PI, "Sensor field of view.");
DEFINE_double(regularizer, 1.0, "Regularization parameter for belief update.");

using namespace radiation;

// Create a globally-defined ExplorerLP
ExplorerLP explorer(FLAGS_num_rows, FLAGS_num_cols, FLAGS_num_sources,
                    FLAGS_regularizer, FLAGS_num_steps, FLAGS_fov,
                    FLAGS_num_samples);

// Create a globally-defined step counter.
unsigned int step_count = 0;

// Initialize OpenGL.
void InitGL() {
  // Set the "clearing" or background color as black/opaque.
  glClearColor(0.0, 0.0, 0.0, 1.0);

  // Set up alpha blending.
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable( GL_BLEND );
}

// Timer callback. Re-render at the specified rate.
void Timer(int value) {
  glutPostRedisplay();
  glutTimerFunc(FLAGS_refresh_rate, Timer, 0);
}

// Reshape the window to maintain the correct aspect ratio.
void Reshape(GLsizei width, GLsizei height) {
  if (height == 0)
    height = 1;

  // Compute aspect ratio fo the new window.
  const GLfloat kAspectRatio =
    static_cast<GLfloat>(width) / static_cast<GLfloat>(height);

  // Set the viewport to cover the new window.
  glViewport(0, 0, width, height);

  // Set the clipping area to be a square in the positive quadrant.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (width >= height) {
    // Larger width than height.
    gluOrtho2D(0.0, static_cast<GLfloat>(FLAGS_num_rows) * kAspectRatio,
               0.0, static_cast<GLfloat>(FLAGS_num_cols));
  } else {
    // Larger height than width.
    gluOrtho2D(0.0, static_cast<GLfloat>(FLAGS_num_rows),
               0.0, static_cast<GLfloat>(FLAGS_num_cols) / kAspectRatio);
  }
}

// Run a single iteration of the exploration algorithm.
void SingleIteration() {
  // Only plan ahead if we haven't exceeded the step count and the entropy
  // is still large.
  if ((FLAGS_iterate_forever || step_count < FLAGS_num_iterations) &&
      (explorer.Entropy() > 1.0)) {
    // Plan ahead.
    std::vector<GridPose2D> trajectory;
    if (!explorer.PlanAhead(trajectory)) {
      std::cout << "error in exploration" << std::endl << std::flush;
      VLOG(1) << "Explorer encountered an error. Skipping this iteration.";
      return;
    }

    // Take a step.
    const double entropy = explorer.TakeStep(trajectory);
    step_count++;
    std::cout << "Entropy after step " << step_count <<
      " is " << entropy << "." << std::endl << std::flush;
  }

  // No matter what, visualize.
  explorer.Visualize();

  // Swap buffers.
  glutSwapBuffers();
}

// Set everything up and go!
int main(int argc, char** argv) {
  // Set up logging.
  google::InitGoogleLogging(argv[0]);

  // Parse flags.
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Set static variables.
  GridPose2D::SetNumRows(FLAGS_num_rows);
  GridPose2D::SetNumCols(FLAGS_num_cols);
  Movement2D::SetAngularStep(FLAGS_angular_step);

  // Set up OpenGL window.
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE);
  glutInitWindowSize(320, 320);
  glutInitWindowPosition(50, 50);
  glutCreateWindow("ExplorerLP");
  glutDisplayFunc(SingleIteration);
  glutReshapeFunc(Reshape);
  glutTimerFunc(0, Timer, 0);
  InitGL();
  glutMainLoop();

  return 0;
}
