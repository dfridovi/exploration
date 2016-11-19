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

#ifndef RADIATION_MOVEMENT_2D_H
#define RADIATION_MOVEMENT_2D_H

#include <vector>
#include <random>

namespace radiation {

class Movement2D {
public:
  ~Movement2D();

  // Default constructor picks a random perturbation dx, dy, da.
  // Alternatively, construct by specifying indices into delta arrays.
  Movement2D();
  Movement2D(unsigned int x_id, unsigned int y_id, unsigned int a_id);

  // Static setters.
  static void SetDeltaXs(const std::vector<double>& delta_xs);
  static void SetDeltaYs(const std::vector<double>& delta_ys);
  static void SetDeltaAngles(const std::vector<double>& delta_as);
  static void SetAngularStep(double angular_step);

  // Getters.
  static unsigned int GetNumDeltaXs();
  static unsigned int GetNumDeltaYs();
  static unsigned int GetNumDeltaAngles();

  unsigned int GetIndexX() const;
  unsigned int GetIndexY() const;
  unsigned int GetIndexAngle() const;

  double GetDeltaX() const;
  double GetDeltaY() const;
  double GetDeltaAngle() const;

private:
  // Static variables. Sets of dx, dy, da, where actually the real change in
  // angle will be angular_step_ * delta_as_. Also a random number generator.
  static std::vector<double> delta_xs_;
  static std::vector<double> delta_ys_;
  static std::vector<double> delta_as_;
  static double angular_step_;

  static std::random_device rd_;
  static std::default_random_engine rng_;

  // Non-static variables. Indices in the static delta vectors.
  unsigned int xx_, yy_, aa_;
}; // struct Movement2D

} // namespace radiation

#endif
