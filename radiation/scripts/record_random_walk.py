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
# Record a dataset of map, successive poses, and corresponding measurements.
# The dataset will contain many such groups.
#
###########################################################################

from grid_pose_2d import GridPose2D
from source_2d import Source2D
from sensor_2d import Sensor2D
from encoding import *

import numpy as np
import math

# File to save to.
maps_file = "maps.csv"
trajectories_file = "trajectories.csv"
measurements_file = "measurements.csv"

# Define hyperparameters.
kNumSimulations = 100000
kNumRows = 5
kNumCols = 5
kNumSources = 3
kNumSteps = 10
kNumAngles = 8
kAngularStep = 2.0 * math.pi / float(kNumAngles)
kSensorParams = {"x" : kNumRows/2,
                 "y" : kNumCols/2,
                 "angle" : 0.0,
                 "fov" : 0.5 * math.pi}
delta_xs = [-1, 0, 1]
delta_ys = [-1, 0, 1]
delta_as = [-kAngularStep, 0, kAngularStep]

# Run the specified number of simulations.
maps = np.zeros(kNumSimulations)
trajectories = np.zeros((kNumSimulations, kNumSteps))
measurements = np.zeros((kNumSimulations, kNumSteps))
for ii in range(kNumSimulations):
    # Generate random sources on the grid.
    sources = []
    for jj in range(kNumSources):
        x = float(np.random.random_integers(0, kNumRows-1)) + 0.5
        y = float(np.random.random_integers(0, kNumCols-1)) + 0.5
        sources.append(Source2D(x, y))

    sensor = Sensor2D(kSensorParams, sources)
    maps[ii] = EncodeMap(kNumRows, kNumCols, sources)

    # Generate a valid trajectory of the given length.
    step_counter = 0
    current_pose = GridPose2D(kNumRows, kNumCols,
                              int(np.random.uniform(0.0, kNumRows)) + 0.5,
                              int(np.random.uniform(0.0, kNumCols)) + 0.5,
                              np.random.uniform(0.0, 2.0 * math.pi))
    while step_counter < kNumSteps:
        dx = np.random.choice(delta_xs)
        dy = np.random.choice(delta_ys)
        da = np.random.choice(delta_as)

        next_pose = GridPose2D.Copy(current_pose)
        if next_pose.MoveBy(dx, dy, da):
            # If a valid move, append to list.
            trajectories[ii, step_counter] = (int(next_pose.x_) +
                                              int(next_pose.y_) * kNumRows +
                                              (int(next_pose.angle_ / kAngularStep)
                                               % kNumAngles) * kNumRows * kNumAngles)
            current_pose = next_pose

            # Get a measurement.
            sensor.ResetPose(current_pose)
            measurements[ii, step_counter] = sensor.Sense()
            step_counter += 1

# Save to disk.
np.savetxt(maps_file, maps, delimiter=",")
np.savetxt(trajectories_file, trajectories, delimiter=",")
np.savetxt(measurements_file, measurements, delimiter=",")
print "Successfully saved to disk."
