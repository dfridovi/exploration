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
from encoding import *

import numpy as np
import math

class Problem:
    def __init__(self, num_rows, num_cols, num_sources, num_steps,
                 angular_step, sensor_params, num_samples, map_prior=None):
        self.num_rows_ = int(num_rows)
        self.num_cols_ = int(num_cols)
        self.num_sources_ = int(num_sources)
        self.num_steps_ = int(num_steps)
        self.angular_step_ = angular_step
        self.sensor_params_ = sensor_params
        self.num_samples_ = int(num_samples)
        self.map_prior_ = map_prior

        # If 'map_prior' is not set, then create a new uniform prior.
        if (self.map_prior_ is None):
            self.map_prior_ = GridMap2D(self.num_rows_, self.num_cols_,
                                        self.num_sources_)

    def GenerateConditionals(self, pose):
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

        # Set the choices for taking steps.
        delta_xs = [-1, 0, 1]
        delta_ys = [-1, 0, 1]
        delta_as = [-self.angular_step_, 0.0, self.angular_step_]

        # Compute the number of possible trajectories, maps, and measurements.
        kNumTrajectories = (len(delta_xs)*len(delta_ys)*len(delta_as))**self.num_steps_
        kNumMaps = (self.num_rows_ * self.num_cols_)**self.num_sources_
        kNumMeasurements = (self.num_sources_ + 1)**self.num_steps_

        # Create empty matrices to store the joint distributions [P_{Z, X}]
        # and [P_{M, Z}] (from which we can estimate what we need).
        zx_joint = np.zeros((kNumMeasurements, kNumTrajectories))
        zm_joint = np.zeros((kNumMeasurements, kNumMaps))

        # Generate a ton of sampled data.
        for ii in range(self.num_samples_):
            # Generate random sources on the grid according to 'map_prior' and
            # compute a corresponding map id number based on which grid cells
            # the sources lie in.
            sources = self.map_prior_.GenerateSources()
            map_id = EncodeMap(self.num_rows_, self.num_cols_, sources)

            # Create a sensor.
            sensor = Sensor2D(self.sensor_params_, sources)

            # Pick a random trajectory starting at the given pose. At each step,
            # get the corresponding measurement.
            current_pose = pose
            delta_sequence = []
            measurements = []

            while len(delta_sequence) < self.num_steps_:
                dx = np.random.choice(delta_xs)
                dy = np.random.choice(delta_ys)
                da = np.random.choice(delta_as)

                next_pose = GridPose2D.Copy(current_pose)
                if next_pose.MoveBy(dx, dy, da):
                    # If a valid move, append to list.
                    delta_sequence.append((dx, dy, da))
                    current_pose = next_pose

                    # Get a measurement.
                    sensor.ResetPose(current_pose)
                    measurements.append(sensor.Sense())

            # Get trajectory and measurement ids.
            trajectory_id = EncodeTrajectory(delta_xs, delta_ys, delta_as,
                                             delta_sequence)
            measurement_id = EncodeMeasurements(self.num_sources_, measurements)

            # Record this sample.
            zx_joint[measurement_id, trajectory_id] += 1.0
            zm_joint[measurement_id, map_id] += 1.0

        # Normalize so that all the rows sum to unity.
        zx_conditional = zx_joint / np.sum(zx_joint, axis=1)[:, None]
        zm_conditional = zm_joint / np.sum(zm_joint, axis=1)[:, None]

        # Compute [h_{M|Z}], the conditional entropy vector.
        def entropy(distribution):
            return sum(map(lambda p : -max(p, 1e-4) * math.log(max(p, 1e-4)),
                           distribution))

        h_conditional = np.asarray(
            map(lambda measurement_id : entropy(zm_conditional[measurement_id, :]),
                range(kNumMeasurements)))

        return (zx_conditional, h_conditional)
