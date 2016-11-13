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
  bool DecodeTrajectory(unsigned int id, unsigned int num_steps,
                        const GridPose2D& initial_pose,
                        std::vector<GridPose2D>& trajectory) {

  }

  // Encode/decode measurements.
  unsigned int EncodeMeasurements(const std::vector<unsigned int>& measurements,
                                  unsigned int max_measurement);
  bool DecodeMeasurements(unsigned int id, unsigned int max_measurement,
                          unsigned int num_measurements,
                          std::vector<unsigned int>& measurements);

  // Encode/decode list of sources (map).
  unsigned int EncodeMap(const std::vector<Source2D>& sources,
                         unsigned int num_rows, unsigned int num_cols);
  bool DecodeMap(unsigned int id, unsigned int num_rows, unsigned int num_cols,
                 unsigned int num_sources, std::vector<Source2D>& sources);

} // namespace radiation

#endif
