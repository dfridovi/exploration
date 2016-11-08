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
from problem import Problem
from encoding import *

from scipy.optimize import linprog
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

class ExplorerLP:
    def __init__(self, num_rows, num_cols, num_sources, num_steps,
                 angular_step, sensor_params, num_samples):
        """ Constructor. """
        self.num_rows_ = num_rows
        self.num_cols_ = num_cols
        self.num_sources_ = num_sources
        self.num_steps_ = num_steps
        self.angular_step_ = angular_step
        self.sensor_params_ = sensor_params
        self.num_samples_ = num_samples

        # Keep track of map, pose, true sources, and past poses.
        self.map_ = GridMap2D(self.num_rows_, self.num_cols_, self.num_sources_)
        self.pose_ = GridPose2D(self.num_rows_, self.num_cols_,
                                0.5 * self.num_rows_, 0.5 * self.num_cols_, 0.0)
        self.sources_ = []
        self.past_poses_ = []

        # Generate random sources on the grid.
        for ii in range(self.num_sources_):
            x = float(np.random.random_integers(0, self.num_rows_-1)) + 0.5
            y = float(np.random.random_integers(0, self.num_cols_-1)) + 0.5
            self.sources_.append(Source2D(x, y))

    def PlanAhead(self):
        """
        Generate the reward vector for the LP formulation of planning and solve.
        """
        problem = Problem(self.num_rows_, self.num_cols_, self.num_sources_,
                          self.num_steps_, self.angular_step_, self.sensor_params_,
                          self.num_samples_, self.map_)

        # Generate LP parameters at the current pose.
        (pzx, hzm, trajectory_ids) = problem.GenerateConditionals(self.pose_)
        pzx = np.asmatrix(pzx)
        hzm = np.asmatrix(hzm).T
        num_trajectories = pzx.shape[1]

        # Set up LP.
        result = linprog(np.asarray((pzx.T * hzm)).ravel(),
                         A_eq=np.ones((1, num_trajectories)),
                         b_eq=np.ones((1, 1)),
                         method='simplex')

        if (not result.success):
            print "Could not find a feasible solution to the LP."
            return []

        # Decode solution into trajectory.
        print "Found information-optimal trajectory."
        trajectory_id = trajectory_ids[np.argmax(result.x)]

        delta_xs = [-1, 0, 1]
        delta_ys = [-1, 0, 1]
        delta_as = [-self.angular_step_, 0.0, self.angular_step_]
        trajectory = DecodeTrajectory(delta_xs, delta_ys, delta_as,
                                      trajectory_id, self.pose_, self.num_steps_)

        return trajectory

    def TakeStep(self, trajectory):
        """ Move one step along this trajectory. """

        # Update list of past poses.
        current_pose = GridPose2D.Copy(self.pose_)
        self.past_poses_.append(current_pose)

        # Update pose.
        self.pose_ = trajectory[0]
        self.sensor_params_["x"] = self.pose_.x_
        self.sensor_params_["y"] = self.pose_.y_
        self.sensor_params_["angle"] = self.pose_.angle_

        # Take scan, and update map.
        sensor = Sensor2D(self.sensor_params_, self.sources_)
        self.map_.Update(sensor)

        # Return entropy.
        return self.map_.Entropy()

    def Visualize(self, title):
        """
        Display a visualization of the current belief state,
        the true locations of the sources, the pose, and the field of view.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")
        ax.set_xticks(np.arange(0, self.map_.belief_.shape[1], 1))
        ax.set_yticks(np.arange(0, self.map_.belief_.shape[0], 1))

        # Plot the belief grid, one square at a time.
        patches = []
        colors = []
        for ii in range(self.map_.belief_.shape[0]):
            for jj in range(self.map_.belief_.shape[1]):
                patch = mpatches.Rectangle((float(ii), float(jj)), 1.0, 1.0)
                patches.append(patch)
                colors.append(self.map_.belief_[ii, jj])

        patch_collection = PatchCollection(patches, cmap=plt.cm.bone, alpha=0.9)
        patch_collection.set_array(np.array(colors))
        ax.add_collection(patch_collection)

        try:
            plt.colorbar(patch_collection)
        except Exception, e:
            pass

        # Overlay the robot position as a circle.
        ax.scatter([self.pose_.x_], [self.pose_.y_],
                    s=[np.pi * 15**2], color="blue", alpha=0.75)

        # Overlay the robot's field of view as a colored wedge.
        fov = self.sensor_params_["fov"]
        upper_bound = self.pose_.angle_ + 0.5 * fov
        lower_bound = self.pose_.angle_ - 0.5 * fov

        wedge = mpatches.Wedge((self.pose_.x_, self.pose_.y_),
                               1.5 * max(self.map_.belief_.shape[0],
                                          self.map_.belief_.shape[1]),
                               180.0 * lower_bound / np.pi,
                               180.0 * upper_bound / np.pi,
                               facecolor="blue",
                               alpha=0.5)
        ax.add_patch(wedge)

        # Overlay all past poses, with their fields of view.
        for ii, pose in enumerate(self.past_poses_):
            past_upper_bound = pose.angle_ + 0.5 * fov
            past_lower_bound = pose.angle_ - 0.5 * fov

            fade = 0.1 * float(ii + 1) / len(self.past_poses_)
            past_wedge = mpatches.Wedge((pose.x_, pose.y_),
                                        1.5 * max(self.map_.belief_.shape[0],
                                                  self.map_.belief_.shape[1]),
                                        180.0 * past_lower_bound / np.pi,
                                        180.0 * past_upper_bound / np.pi,
                                        facecolor="green",
                                        alpha=fade)
            ax.add_patch(past_wedge)

            ax.scatter([pose.x_], [pose.y_],
                       s=[np.pi * 15**2], color="green", alpha=fade)

        # Overlay the position of all sources.
        for source in self.sources_:
            ax.scatter([source.x_], [source.y_],
                        s=[np.pi * 7.5**2], color="red", alpha=0.75)


        plt.title(title)
        ax.set_xlim([-0.5, self.map_.belief_.shape[0] + 0.5])
        ax.set_ylim([-0.5, self.map_.belief_.shape[1] + 0.5])
        ax.set_xticks(range(self.map_.belief_.shape[0] + 1))
        ax.set_yticks(range(self.map_.belief_.shape[1] + 1))

        plt.show()
