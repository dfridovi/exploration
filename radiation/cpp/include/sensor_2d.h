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
// Defines a 2D radiation sensor, which can sense sources in its fied of view.
// Returns a noiseless count of the sources in view.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef RADIATION_SENSOR_2D_H
#define RADIATION_SENSOR_2D_H

#include <source_2d.h>
#include <grid_pose_2d.h>

#include <vector>

namespace radiation {

class Sensor2D {
public:
  Sensor2D(double x, double y, double a, double fov);
  Sensor2D(const GridPose2D& pose, double fov);
  ~Sensor2D();

  // Getters.
  double GetX() const;
  double GetY() const;
  double GetAngle() const;

  unsigned int GetIndexX() const;
  unsigned int GetIndexY() const;

  // Move to the given location and orientation.
  void MoveTo(double x, double y, double a);

  // Sense the specified sources.
  unsigned int Sense(const std::vector<Source2D>& sources) const;

  // Check if a source or voxel is in view.
  bool SourceInView(const Source2D& source) const;
  bool VoxelInView(unsigned int ii, unsigned int jj) const;

private:
  // Position and orientation angle.
  double x_, y_, a_;

  // Field of view.
  const double fov_;

}; // struct Source2D

} // namespace radiation

#endif
