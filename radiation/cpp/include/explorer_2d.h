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
// Exploration on a 2D grid. All explorers will inherit from this class, which
// provides basic functionality like taking steps along trajectories,
// computing entropy, and visualization.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIATION_EXPLORER_2D_H
#define RADIATION_EXPLORER_2D_H

#include <source_2d.h>
#include <sensor_2d.h>
#include <grid_map_2d.h>
#include <grid_pose_2d.h>
#include <movement_2d.h>
#include <encoding.h>

#include <Eigen/Core>
#include <vector>

namespace radiation {

class Explorer2D {
 public:
  Explorer2D(unsigned int num_rows, unsigned int num_cols,
             unsigned int num_sources, double regularizer,
             double fov, const std::vector<Source2D>& sources,
             const GridPose2D& initial_pose);
  Explorer2D(unsigned int num_rows, unsigned int num_cols,
             unsigned int num_sources, double regularizer,
             double fov);
  virtual ~Explorer2D();

  // Plan a new trajectory. Must implement this in all derived classes.
  virtual bool PlanAhead(std::vector<GridPose2D>& trajectory) = 0;

  // Take a step along the given trajectory. Return resulting entropy.
  double TakeStep(const std::vector<GridPose2D>& trajectory);

  // Compute map entropy.
  double Entropy() const;

  // Visualize the current belief state.
  void Visualize() const;

 protected:
  // Problem parameters.
  double fov_;

  // Map, pose, and sources.
  GridMap2D map_;
  GridPose2D pose_;
  std::vector<Source2D> sources_;

  // List of past poses.
  std::vector<GridPose2D> past_poses_;
}; // class Explorer2D

} // namespace radiation

#endif
