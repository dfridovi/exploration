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
# Defines a grid map in 2D, with the specified dimensions.
# There is no true length scale in this structure; angles however are
# assumed to be metricly accurate.
#
# Implicitly, the values in the grid define a 2D non-homogeneous Poisson
# Process, which allows the random generation of possible sources. This
# is used to simulate measurements from arbitrary locations based on
# the current belief state.
#
###########################################################################

import numpy as np
from numpy import linalg as LA
import math

class GridMap2D:
    def __init__(self, nrows, ncols, k):
        """
        Constructor. Takes in dimensions and number of sources k
        and generates a uniform prior.
        """
        self.belief_ = np.ones((nrows, ncols), dtype=np.float) * k / (nrows*ncols)
        self.k_ = k

    def Update(self, sensor):
        """
        Update belief about the world, given a sensor (with associated
        paramters, including position and orientation).

        For all voxels in range, create a rate update which is uniform
        and sums to the measurement value, then compute a weighted average
        at each point, trusting lower values more than higher ones.
        """
        measurement = sensor.Sense()
        if measurement > self.k_:
            print "Measured too many sources. Did not update."
            return False

        # Identify voxels that are in range.
        update = np.copy(self.belief_)
        in_view_mask = np.zeros((self.belief_.shape), dtype=np.bool)
        in_view_count = 0
        for ii in range(self.belief_.shape[0]):
            for jj in range(self.belief_.shape[1]):
                if sensor.VoxelInView(ii, jj):
                    in_view_mask[ii, jj] = True
                    in_view_count += 1

        if (in_view_count) == 0:
            assert measurement == 0
            return True

        # Set update array.
        update[np.where(in_view_mask)] = float(measurement) / in_view_count

        # Catch situation where we might divide by zero.
        if in_view_count < update.shape[0] * update.shape[1]:
            update[np.where(np.invert(in_view_mask))] = (float(self.k_ - measurement) /
                                        (update.shape[0]*update.shape[1] - in_view_count))

        # Perform update.
        update_weights_top = np.abs(self.belief_)
        update_weights_bottom = np.abs(update) + update_weights_top

        # Only update where belief is sufficiently large.
        valid_indices = np.where(update_weights_top > 1e-8)

        update_weights = np.zeros(update.shape)
        update_weights[valid_indices] = np.divide(update_weights_top[valid_indices],
                                                  update_weights_bottom[valid_indices])
        self.belief_ = ((1.0 - update_weights) * self.belief_ +
                        np.multiply(update_weights, update))
        return True

    def SimulateTrajectory(self, sensor, trajectory, niters=1):
        """
        Return expected map entropy after taking scans at each GridPose2D
        in the trajectory.

        Expectation is based on Monte Carlo simulation using the specified
        number of iterations.
        """
        entropy_total = 0

        # Save the current belief state.
        current_belief = np.copy(self.belief_)

        # Monte Carlo simulation.
        for ii in range(niters):
            for pose in trajectory:
                # Set sensor pose.
                sensor.ResetPose(pose)

                # Simulate sensor measurement.
                assert self.Update(sensor)

            # Once all poses have been simulated, compute entropy.
            entropy_total += self.Entropy()

            # Reset belief state.
            self.belief_ = np.copy(current_belief)

        # Divide out by number of iterations.
        return entropy_total / float(niters)

    def SimulatePose(self, sensor, niters=1):
        """
        Return expected map entropy after receiving a measurement from
        the specified location/orientation. Expectation is based on
        Monte Carlo simulation using the specified number of iterations.
        """
        entropy_total = 0

        # Save the current belief state.
        current_belief = np.copy(self.belief_)

        # Iterate the specified number of times, accumulating total entropy.
        for ii in range(niters):
            assert self.Update(sensor)
            entropy_total += self.Entropy()

            # Reset belief state.
            self.belief_ = np.copy(current_belief)

        # Divide out by number of iterations.
        return entropy_total / float(niters)

    def Entropy(self):
        """
        Compute the entropy of the map. Since we model each voxel as a
        Bernoulli variable, independent of all others, we compute joint entropy
        as the sum of binary entropies at each voxels.
        """
        total_entropy = 0
        for ii in range(self.belief_.shape[0]):
            for jj in range(self.belief_.shape[1]):
                total_entropy += self.BernoulliEntropy(self.belief_[ii, jj])

        return total_entropy

    def BernoulliEntropy(self, p):
        """
        Binary entropy function (in nats) of a Bernoulli variable with
        parameter p.
        """
        assert (p >= 0.0) and (p <= 1.0)

        # For numerical stability, return 0 if p is nearly 0 or 1.
        if (p < 1e-8) or (p > 1.0 - 1e-8):
            return 0.0

        # Otherwise, use the definition of entropy.
        return -p * math.log(p) - (1.0 - p) * math.log(1.0 - p)

    def PoissonEntropy(self, rate):
        """
        Approximate entropy (in nats) of a Poisson variable with the
        given rate.
        """
        if rate < 1e-8:
            return 0.0

        sum_max = min(6, 2*int(round(rate)))
        entropy = rate * (1.0 - math.log(rate))

        # Add on sum_max extra terms.
        extra_terms = 0
        rate_ii = rate
        fact_ii = 1
        for ii in range(2, sum_max + 1):
            rate_ii *= rate
            fact_ii *= float(ii)
            extra_terms += rate_ii * math.log(fact_ii) / fact_ii

        # Scale extra terms and add to entropy.
        entropy += math.exp(-rate) * extra_terms
        return entropy
