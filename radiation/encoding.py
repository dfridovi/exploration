"""
Copyright (c) 2015, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

   3. Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

Please contact the author(s) of this library if you have any questions.
Authors: David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
"""

###########################################################################
#
# Encoders and decoders for measurements, trajectories, and maps.
#
###########################################################################

from grid_pose_2d import GridPose2D
from grid_map_2d import GridMap2D
from source_2d import Source2D
from sensor_2d import Sensor2D

import numpy as np
import math

# Encode a list of (dx, dy, da) tuples in an integer from 0 to
# (len(delta_xs) * len(delta_ys) * len(delta_as))**len(delta_sequence) - 1.
def EncodeTrajectory(delta_xs, delta_ys, delta_as, delta_sequence):
    base = len(delta_xs) * len(delta_ys) * len(delta_as)
    trajectory_id = 0

    for ii, delta in enumerate(delta_sequence):
        x_id = delta_xs.index(delta[0])
        y_id = delta_ys.index(delta[1])
        a_id = delta_as.index(delta[2])

        delta_id = x_id + y_id*len(delta_xs) + a_id*len(delta_xs)*len(delta_ys)
        trajectory_id += delta_id * base**ii

    return trajectory_id

# Decode a trajectory id into a list of GridPose2Ds.
def DecodeTrajectory(delta_xs, delta_ys, delta_as, trajectory_id,
                     initial_pose, num_steps):
    base = len(delta_xs) * len(delta_ys) * len(delta_as)
    trajectory = []

    current_pose = initial_pose
    while trajectory_id > 0:
        remainder = trajectory_id  % base

        # Convert remainder to delta tuple (dx, dy, da).
        x_id = remainder % len(delta_xs)
        y_id = ((remainder - x_id) / len(delta_xs)) % len(delta_ys)
        a_id = ((remainder - x_id - y_id*len(delta_xs)) /
                (len(delta_xs)*len(delta_ys)))

        dx = delta_xs[x_id]
        dy = delta_ys[y_id]
        da = delta_as[a_id]

        # Append to 'trajectory'.
        next_pose = GridPose2D.Copy(current_pose)
        assert next_pose.MoveBy(dx, dy, da)
        trajectory.append(next_pose)

        # Update 'trajectory_id'.
        trajectory_id = (trajectory_id - remainder) / base

        # Reset 'current_pose'.
        current_pose = next_pose

    # If not the right length, that means that the last remainders were 0.
    # Update 'trajectory' accordingly.
    while len(trajectory) < num_steps:
        next_pose = GridPose2D.Copy(current_pose)
        assert next_pose.MoveBy(delta_xs[0], delta_ys[0], delta_as[0])

        trajectory.append(next_pose)
        current_pose = next_pose

    return trajectory

# Encode a list of measurements in an integer.
def EncodeMeasurements(max_measurement, measurement_sequence):
    base = max_measurement + 1
    measurement_id = 0

    for ii, measurement in enumerate(measurement_sequence):
        measurement_id += measurement * base**ii

    return measurement_id

# Decode a measurement id into a list of measurements.
def DecodeMeasurements(max_measurement, measurement_id, num_measurements):
    base = max_measurement + 1
    measurements = []

    while measurement_id > 0:
        remainder = measurement_id % base
        measurements.append(remainder)

        # Update 'measurement_id'.
        measurement_id = (measurement_id - remainder) / base

    # If not enough measurements, the rest must have been zero.
    while len(measurements) < num_measurements:
        measurements.append(0)

    return measurements

# Encode a map (list of sources) as an integer.
def EncodeMap(num_rows, num_cols, sources):
    base = num_rows * num_cols
    map_id = 0

    for ii, source in enumerate(sources):
        source_id = int(source.x_) + int(source.y_) * num_rows
        map_id += source_id * base**ii

    return map_id

# Decode a map id into a list of sources.
def DecodeMap(num_rows, num_cols, map_id, num_sources):
    base = num_rows * num_cols
    sources = []

    while map_id > 0:
        remainder = map_id % base

        # Unpack remainder into (x, y) coordinates of a source.
        source = Source2D(float(remainder % num_rows) + 0.5,
                          float(remainder // num_rows) + 0.5)
        sources.append(source)

        # Update 'map_id'.
        map_id = (map_id - remainder) / base

    # If not enough sources, the remainders must have been zero.
    while len(sources) < num_sources:
        sources.append(Source2D(0.5, 0.5))

    return sources
