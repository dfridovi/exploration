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
# Defines a problem setup, parameterized by grid size, number of sources,
# number of steps, angular step size, sensor parameters, and number of
# samples to estimate distributions from.
#
# The key functionality here is to generate the conditional distribution
# matrix [P_{Z|X}] and the conditional entropy vector [h_{M|Z}].
#
###########################################################################

from grid_pose_2d import GridPose2D
from grid_map_2d import GridMap2D
from source_2d import Source2D
from sensor_2d import Sensor2D

import numpy as np

class Problem:
    def __init__(self, num_rows, num_cols, num_sources,
                 num_steps, angular_step, sensor_params, num_samples):
        self.num_rows_ = int(num_rows)
        self.num_cols_ = int(num_cols)
        self.num_sources_ = int(num_sources)
        self.num_steps_ = int(num_steps)
        self.angular_step_ = angular_step
        self.sensor_params_ = sensor_params
        self.num_samples_ = int(num_samples)

    def GenerateConditionalDistribution(self, pose):
        """
        Generate conditional distribution matrix [P_{Z|X}] and conditional
        entropy vector [h_{M|Z}] by starting from the specified pose and
        generating a bunch of random legal trajectories, and for each one
        generating a random map and a measurement.

        The (i,j)-entry of [P_{Z|X}] is the normalized frequency of observing
        measurement i given that the trajectory chosen was j. This is computed
        by first estimating the joint distribution of Z and X and then
        normalizing each row so that it sums to unity.

        The i-entry of [h_{M|Z}] is the entropy of M given that we observed
        measurement i, starting from the given pose. This is computed by first
        estimating the joint distribution of M and Z, and then looking at the
        row where Z = i.
        """

        # Generate a ton of sampled data.
        samples = []
        for ii in range(self.num_samples_):

            # Pick a random trajectory from the given pose.
            current_pose = pose
            trajectory = []

            while len(trajectory) < self.num_steps_:
                delta_x = np.random.random_integers(-1, 1)
                delta_y = np.random.random_integers(-1, 1)
                delta_angle = (self.angular_step_ *
                               float(np.random.random_integers(-1, 1)))

                next_pose = GridPose2D.Copy(current_pose)
                if next_pose.MoveBy(delta_x, delta_y, delta_angle):
                    trajectory.append(next_pose)
                    current_pose = next_pose

            # Generate random sources on the grid and compute a corresponding
            # map id number based on which grid cells the sources lie in.
            sources = []
            map_id = 0
            for ii in range(k):
                x = math.floor(
                    np.random.uniform(0.0, float(self.num_rows_))) + 0.5
                y = math.floor(
                    np.random.uniform(0.0, float(self.num_cols_))) + 0.5

                sources.append(Source2D(x, y))
                map_id += ((self.num_cols_ * int(x) + int(y)) *
                           (self.num_rows_ * self.num_cols_)**ii)

            # Create a sensor.
            sensor = Sensor2D(self.sensor_params_, sources)

            # Walk the trajectory and obtain measurements. Compute a measurement
            # id number based on which measurement was obtained at which step.
            measurement_id = 0
            for ii, step in enumerate(trajectory):
                sensor.ResetPose(step)
                measurement = sensor.Sense()

                # Update measurement id.
                measurement_id += measurement * self.num_sources_**ii

            # Record this sample.
            samples.append((trajectory_id, map_id, measurement_id))
