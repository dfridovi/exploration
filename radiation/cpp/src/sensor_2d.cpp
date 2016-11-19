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

#include <sensor_2d.h>

#include <math.h>

namespace radiation {

  Sensor2D::~Sensor2D() {}
  Sensor2D::Sensor2D(double x, double y, double a, double fov)
    : x_(x), y_(y), a_(a), fov_(fov) {}
  Sensor2D::Sensor2D(const GridPose2D& pose, double fov)
    : x_(pose.GetX()), y_(pose.GetY()), a_(pose.GetAngle()), fov_(fov) {}

  // Getters.
  double Sensor2D::GetX() const { return x_; }
  double Sensor2D::GetY() const { return y_; }
  double Sensor2D::GetAngle() const { return a_; }

  unsigned int Sensor2D::GetIndexX() const {
    return static_cast<unsigned int>(x_);
  }

  unsigned int Sensor2D::GetIndexY() const {
    return static_cast<unsigned int>(y_);
  }

  // Move to the given location and orientation.
  void Sensor2D::MoveTo(double x, double y, double a) {
    x_ = x;
    y_ = y;
    a_ = a;
  }

  // Sense the specified sources. Count the number in view.
  unsigned int Sensor2D::Sense(const std::vector<Source2D>& sources) const {
    unsigned int count = 0;

    for (const auto& source : sources) {
      if (SourceInView(source))
        count++;
    }

    return count;
  }

  // Check if a source or voxel is in view.
  bool Sensor2D::VoxelInView(unsigned int ii, unsigned int jj) const {
    return SourceInView(Source2D(ii, jj));
  }

  bool Sensor2D::SourceInView(const Source2D& source) const {
    // Get unit vector to source.
    double dx = source.GetX() - x_;
    double dy = source.GetY() - y_;
    double norm = sqrt(dx * dx + dy * dy);

    if (norm < 1e-8)
      return true;

    dx /= norm;
    dy /= norm;

    // In view if the angle of this vector int he plane is within half the
    // field of view of our current angle.
    double angle_to_source = acos(dx * cos(a_) + dy * sin(a_));
    if (angle_to_source < 0.5 * fov_)
      return true;

    return false;
  }

} // namespace radiation
