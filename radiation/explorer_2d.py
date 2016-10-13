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
# Exploration in 2D on a grid. Tries to find k the radiation sources in
# as few steps as possible by taking the step which maximizes mutual
# information between expected measurements and the true map.
#
###########################################################################

from grid_pose_2d import GridPose2D
from grid_map_2d import GridMap2D
from source_2d import Source2D
from sensor_2d import Sensor2D

import numpy as np
import matplotlib.pyplot as plt

class Explorer2D:
    def __init__(self, nrows, ncols, k, angular_step, sensor_params):
        """
        Constructor. Takes in dimensions, number of sources, resolution,
        and sensor parameters.
        """
        self.angular_step_ = angular_step
        self.sensor_params_ = sensor_params
        self.map_ = GridMap2D(nrows, ncols, k)
        self.pose_ = GridPose2D(nrows, ncols, 0.5 * nrows, 0.5 * ncols, 0.0)
        self.sources_ = []

        for ii in range(k):
            x = np.random.uniform(0.0, float(nrows))
            y = np.random.uniform(0.0, float(ncols))
            self.sources_.append(Source2D(x, y))

    def PlanAhead(self, nsteps, ntrajectories, niters):
        """
        Simulate a bunch of random trajectories and return the best one.
        """
        best_trajectory = []
        best_entropy = float("inf")

        for ii in range(ntrajectories):
            current_pose = GridPose2D.Copy(self.pose_)
            trajectory = []

            # Choose a random trajectory.
            while len(trajectory) < nsteps:
                delta_x = np.random.random_integers(-1, 1)
                delta_y = np.random.random_integers(-1, 1)
                delta_angle = (self.angular_step_ *
                               float(np.random.random_integers(-1, 1)))

                if current_pose.MoveBy(delta_x, delta_y, delta_angle):
                    trajectory.append(current_pose)

            # Compute entropy.
            sensor = Sensor2D(self.sensor_params_, self.sources_)
            entropy = self.map_.SimulateTrajectory(sensor, trajectory, niters)

            # Compare to best.
            if entropy < best_entropy:
                best_trajectory = trajectory
                best_entropy = entropy

        # Return best.
        return best_trajectory

    def TakeStep(self, trajectory):
        """ Move one step along this trajectory. """

        # Update pose.
        self.pose_ = trajectory[0]

        # Take scan, and update map.
        self.sensor_params_["x"] = self.pose_.x_
        self.sensor_params_["y"] = self.pose_.y_
        self.sensor_params_["angle"] = self.pose_.angle_
        sensor = Sensor2D(self.sensor_params_, self.sources_)
        self.map_.Update(sensor)

        # Return entropy.
        return self.map_.Entropy()

    def Visualize(self):
        """
        Display a visualization of the current belief state,
        the true locations of the sources, the pose, and the field of view.
        """
        plt.imshow(self.map_.belief_, cmap=plt.cm.bone)
        plt.colorbar()
        plt.show()
