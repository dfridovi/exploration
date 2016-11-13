
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
// Defines a pose on the 2D grid.
//
///////////////////////////////////////////////////////////////////////////////

#include "grid_pose_2d.h"

#include <math.h>

namespace radiation {

  GridPose2D::~GridPose2D() {}
  GridPose2D::GridPose2D(double x, double y, double a)
    : x_(x), y_(y), a_(a) {}

  // Static setters.
  void GridPose2D::SetNumRows(unsigned int num_rows) { num_rows_ = num_rows; }
  void GridPose2D::SetNumCols(unsigned int num_cols) { num_cols_ = num_cols; }

  // Getters.
  double GridPose2D::GetX() const { return x_; }
  double GridPose2D::GetY() const { return y_; }
  double GridPose2D::GetAngle() const { return a_; }

  unsigned int GridPose2D::GetX() const { return static_cast<unsigned int>(x_); }
  unsigned int GridPose2D::GetY() const { return static_cast<unsigned int>(y_); }

  // Move by the given amount if it is legal.
  bool GridPose2D::MoveBy(const Movement2D& movement) {
    double new_x = x_ + movement.GetDeltaX();
    double new_y = y_ + movement.GetDeltaY();
    double new_a = a_ + movement.GetDeltaAngle();

    // Catch going out of bounds.
    if ((new_x < 0.0) || (new_x > num_rows_) ||
        (new_y < 0.0) || (new_y > num_cols_))
      return false;

    // Not going out of bounds, so update coordinates.
    x_ = new_x;
    y_ = new_y;
    a_ = new_a;

    // Make sure angle is in [0, 2 pi).
    if (a_ < 0.0)
      a_ += M_PI + M_PI;
    else if (a_ > M_PI + M_PI)
      a_ -= M_PI + M_PI;
    return true;
  }

} // namespace radiation
