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
// Exploration on a 2D grid. Tries to find the specified number of radiation
// sources (located at random lattice points) by following a random walk.
//
///////////////////////////////////////////////////////////////////////////////

#include <explorer_rw.h>

#include <glog/logging.h>
#include <random>
#include <string>
#include <math.h>

namespace radiation {

// Constructor/destructor.
ExplorerRW::~ExplorerRW() {}
ExplorerRW::ExplorerRW(unsigned int num_rows, unsigned int num_cols,
                       unsigned int num_sources, double regularizer,
                       double fov, const std::vector<Source2D>& sources,
                       const GridPose2D& initial_pose)
  : Explorer2D(num_rows, num_cols, num_sources, regularizer, fov,
               sources, initial_pose) {}
ExplorerRW::ExplorerRW(unsigned int num_rows, unsigned int num_cols,
                       unsigned int num_sources, double regularizer,
                       double fov)
  : Explorer2D(num_rows, num_cols, num_sources, regularizer, fov) {}

// Plan a new trajectory. Take a completely random (legal) step.
bool ExplorerRW::PlanAhead(std::vector<GridPose2D>& trajectory) {
  trajectory.clear();

  // Pick random movements until one is legal.
  GridPose2D next_pose = pose_;
  while (!next_pose.MoveBy(Movement2D()));

  // Store this next pose.
  trajectory.push_back(next_pose);

  return true;
}

} // namespace radiation
