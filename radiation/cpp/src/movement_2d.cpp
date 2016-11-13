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
// Defines a movement on a 2D grid.
//
///////////////////////////////////////////////////////////////////////////////

#include "movement_2d.h"

#include <glog/logging.h>

namespace radiation {

  Movement2D::~Movement2D() {}
  Movement2D::Movement2D() {
    // Choose each index uniformly from the appropriate set.
    std::uniform_int_distribution<unsigned int> unif_x(0, delta_xs_.size() - 1);
    xx_ = unif_x(rng_);

    std::uniform_int_distribution<unsigned int> unif_y(0, delta_ys_.size() - 1);
    yy_ = unif_y(rng_);

    std::uniform_int_distribution<unsigned int> unif_a(0, delta_as_.size() - 1);
    aa_ = unif_a(rng_);
  }
  Movement2D::Movement2D(unsigned int x_id, unsigned int y_id, unsigned int a_id)
    : xx_(x_id), yy_(y_id), aa_(a_id) {
    CHECK((xx_ >= 0) && (xx_ < delta_xs_.size()));
    CHECK((yy_ >= 0) && (yy_ < delta_ys_.size()));
    CHECK((aa_ >= 0) && (aa_ < delta_as_.size()));
  }

  // Static setters.
  void Movement2D::SetDeltaXs(const std::vector<double>& delta_xs) {
    delta_xs_.clear();
    for (const auto& dx : delta_xs)
      delta_xs_.push_back(dx);
  }

  void Movement2D::SetDeltaYs(const std::vector<double>& delta_ys) {
    delta_ys_.clear();
    for (const auto& dy : delta_ys)
      delta_ys_.push_back(dy);
  }

  void Movement2D::SetDeltaAngles(const std::vector<double>& delta_as) {
    delta_as_.clear();
    for (const auto& da : delta_as)
      delta_as_.push_back(da);
  }

  // Getters.
  unsigned int Movement2D::GetNumDeltaXs() const { return delta_xs_.size(); }
  unsigned int Movement2D::GetNumDeltaYs() const { return delta_ys_.size(); }
  unsigned int Movement2D::GetNumDeltaAngles() const {
    return delta_as_.size();
  }

  unsigned int Movement2D::GetIndexX() const { return xx_; }
  unsigned int Movement2D::GetIndexY() const { return yy_; }
  unsigned int Movement2D::GetIndexAngle() const { return aa_; }

  double Movement2D::GetDeltaX() const { return delta_xs_[xx_]; }
  double Movement2D::GetDeltaY() const { return delta_ys_[yy_]; }
  double Movement2D::GetDeltaAngle() const {
    return angular_step_ * delta_as_[aa_];
  }

} // namespace radiation
