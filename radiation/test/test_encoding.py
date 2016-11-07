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
# Unit tests for encoding/decoding functions.
#
###########################################################################

import numpy as np
import math

from radiation.source_2d import Source2D
from radiation.encoding import *

""" Test encoding/decoding for sources. """
def test_map():
    kNumRows = 10
    kNumCols = 10
    kNumSources = 5

    # Generate a bunch of random sources.
    sources = []
    for ii in range(kNumSources):
        sources.append(Source2D(np.random.random_integers(0, kNumRows-1) + 0.5,
                                np.random.random_integers(0, kNumCols-1) + 0.5))

    # Encode these sources.
    map_id = EncodeMap(kNumRows, kNumCols, sources)

    # Decode the id back into a list of sources.
    decoded_sources = DecodeMap(kNumRows, kNumCols, map_id, kNumSources)

    # Check that each source matches.
    assert len(sources) == len(decoded_sources)
    for original, decoded in zip(sources, decoded_sources):
        dx = original.x_ - decoded.x_
        dy = original.y_ - decoded.y_
        assert math.sqrt(dx*dx + dy*dy) < 1e-8

""" Test encoding/decoding for measurements."""
def test_measurements():
    kMaxMeasurement = 10
    kNumMeasurements = 5

    # Generate a bunch of random measurements.
    measurements = []
    for ii in range(kNumMeasurements):
        measurements.append(np.random.random_integers(0, kMaxMeasurement))

    # Encode these measurments.
    measurement_id = EncodeMeasurements(kMaxMeasurement, measurements)

    # Decode the id back into a list of measurements.
    decoded_measurements = DecodeMeasurements(kMaxMeasurement, measurement_id,
                                              kNumMeasurements)

    # Check that each measurement matches.
    assert len(measurements) == len(decoded_measurements)
    for original, decoded in zip(measurements, decoded_measurements):
        assert original == decoded

""" Test encoding/decoding for trajectories. """
def test_trajectories():
    kNumSteps = 5
    kNumRows = 10
    kNumCols = 10
    delta_xs = [-1, 0, 1]
    delta_ys = [-1, 0, 1]
    delta_as = [-0.125 * math.pi, 0.0, 0.125 * math.pi]

    # Create a random initial pose.
    initial_pose = GridPose2D(kNumRows, kNumCols,
                              float(np.random.random_integers(0, kNumRows-1)) + 0.5,
                              float(np.random.random_integers(0, kNumCols-1)) + 0.5,
                              np.random.uniform(0.0, 2.0 * math.pi))

    # Generate a sequence of movements and the associated trajectory.
    delta_sequence = []
    trajectory = []

    current_pose = initial_pose
    while len(trajectory) < kNumSteps:
        dx = np.random.choice(delta_xs)
        dy = np.random.choice(delta_ys)
        da = np.random.choice(delta_as)

        next_pose = GridPose2D.Copy(current_pose)
        if next_pose.MoveBy(dx, dy, da):
            delta_sequence.append((dx, dy, da))
            trajectory.append(next_pose)
            current_pose = next_pose

    # Encode 'delta_sequence' in an integer.
    trajectory_id = EncodeTrajectory(delta_xs, delta_ys, delta_as, delta_sequence)

    # Decode id back to a trajectory.
    decoded_trajectory = DecodeTrajectory(delta_xs, delta_ys, delta_as, trajectory_id,
                                          initial_pose, kNumSteps)

    # Check that the decoded and actual trajectories match.
    assert len(trajectory) == len(decoded_trajectory)
    for original, decoded in zip(trajectory, decoded_trajectory):
        dx = original.x_ - decoded.x_
        dy = original.y_ - decoded.y_
        da = original.angle_ - decoded.angle_

        assert math.sqrt(dx*dx + dy*dy) < 1e-8
        assert abs(da) < 1e-8
