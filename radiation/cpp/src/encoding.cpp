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
// Encode and decode measurements, trajectories, and maps (list of sources).
//
///////////////////////////////////////////////////////////////////////////////

#include "encoding.h"

namespace radiation {

  // Encode a sequence of movements as an unsigned integer.
  unsigned int EncodeTrajectory(const std::vector<Movement2D>& movements) {
    const unsigned int base =
      Movement2D::GetNumDeltaXs() * Movement2D::GetNumDeltaYs() *
      Movement2D::GetNumDeltaAngles();

    unsigned int id = 0;
    unsigned int place_value = 1;

    for (size_t ii = 0; ii < movements.size(); ii++) {
      const unsigned int x_id = movements[ii].GetIndexX();
      const unsigned int y_id = movements[ii].GetIndexY();
      const unsigned int a_id = movements[ii].GetIndexAngle();

      // Compute a unique identifier for this movement.
      const unsigned int step_id =
        x_id + y_id * Movement2D::GetNumDeltaXs() +
        a_id * Movement2D::GetNumDeltaXs() * Movement2D::GetNumDeltaYs();

      // Update the overall trajectory id, and the place value.
      id += step_id * place_value;
      place_value *= base;
    }

    return id;
  }

  // Decode a trajectory id into a sequence of poses.
  void DecodeTrajectory(unsigned int id, unsigned int num_steps,
                        const GridPose2D& initial_pose,
                        std::vector<GridPose2D>& trajectory) {
    trajectory.clear();
    const unsigned int base =
      Movement2D::GetNumDeltaXs() * Movement2D::GetNumDeltaYs() *
      Movement2D::GetNumDeltaAngles();

    GridPose2D current_pose = initial_pose;
    while (id > 0) {
      remainder = id % base;

      // Convert remainder to Movement2D.
      const unsigned int x_id = remainder % Movement2D::GetNumDeltaXs();
      const unsigned int y_id =
        (remainder / Movement2D::GetNumDeltaXs()) % Movement2D::GetNumDeltaYs();
      const unsigned int a_id =
        remainder / (Movement2D::GetNumDeltaXs() * Movement2D::GetNumDeltaYs());

      const Movement2D step(x_id, y_id, a_id);

      // Append to trajectory.
      GridPose2D next_pose = current_pose;
      CHECK(next_pose.MoveBy(step));
      trajectory.push_back(next_pose);

      // Update 'id'.
      id /= base;

      // Reset 'current_pose'.
      current_pose = next_pose;
    }

    // If the above does not recover enough poses, the last remainders were zero.
    // Update 'trajectory' accordingly.
    while (trajectory.size() < num_steps) {
      GridPose2D next_pose = current_pose;
      CHECK(next_pose.MoveBy(Movement2D(0, 0, 0)));

      trajectory.push_back(next_pose);
      current_pose = next_pose;
    }

    return trajectory;
  }

  // Encode a sequence of measurements in an unsigned integer.
  unsigned int EncodeMeasurements(const std::vector<unsigned int>& measurements,
                                  unsigned int max_measurement) {
    const unsigned int base = max_measurement + 1;

    unsigned int id = 0;
    unsigned int place_value = 1;

    for (size_t ii = 0; ii < measurements.size(); ii++) {
      id += measurements[ii] * place_value;
      place_value *= base;
    }

    return id;
  }

  // Decode a measurement id into a list of measurements.
  void DecodeMeasurements(unsigned int id, unsigned int max_measurement,
                          unsigned int num_measurements,
                          std::vector<unsigned int>& measurements) {
    measurements.clear();
    const unsigned int base = max_measurement + 1;

    while (id > 0) {
      const unsigned int remainder = id % base;
      measurements.push_back(remainder);

      // Update 'id'.
      id /= base;
    }

    // If not enough measurements, the rest must have been zero.
    while (measurements.size() < num_measurements) {
      measurements.push_back(0);
    }
  }

  // Encode a list of sources (map) as an unsigned integer.
  unsigned int EncodeMap(const std::vector<Source2D>& sources,
                         unsigned int num_rows, unsigned int num_cols) {
    const unsigned int base = num_rows * num_cols;

    unsigned int id = 0;
    unsigned int place_value = 1;

    for (size_t ii = 0; ii < sources.size(); ii++) {
      const unsigned int source_id =
        sources[ii].GetIndexX() + sources[ii].GetIndexY() * num_rows;

      id += source_id * place_value;
      place_value *= base;
    }

    return id;
  }

  // Decode a map id into a list of sources.
  void DecodeMap(unsigned int id, unsigned int num_rows, unsigned int num_cols,
                 unsigned int num_sources, std::vector<Source2D>& sources) {
    sources.clear();
    const unsigned int base = num_rows * num_cols;

    while (id > 0) {
      const unsigned int remainder = id % base;

      // Unpack remainder into a source.
      const unsigned int x_id = remainder % num_rows;
      const unsigned int y_id = remainder / num_rows;
      sources.push_back(Source2D(x_id, y_id));

      // Update 'id'.
      id /= base;
    }

    // If not enough sources, the remainders must have been zero.
    while (sources.size() < num_sources) {
      sources.push_back(Source2D(0, 0));
    }
  }

} // namespace radiation
