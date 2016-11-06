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

from source_2d import Source2D

import numpy as np
from numpy import linalg as LA
import math

import copy
from sets import ImmutableSet
from scipy.optimize import least_squares

class GridMap2D:
    def __init__(self, nrows, ncols, k, alpha=1.0):
        """
        Constructor. Takes in dimensions and number of sources k
        and generates a uniform prior. Alpha is a regularization
        parameter used in least squares optimization.
        """
        self.belief_ = np.ones((nrows, ncols), dtype=np.float) * k / (nrows*ncols)
        self.k_ = k
        self.alpha_ = alpha

        # Also, keep track of which voxels have been viewed at each point
        # in time, and the associated measurement.
        self.viewed_lists_ = []
        self.viewed_sets_ = []
        self.measurements_ = []

    def GenerateSources(self):
        """ Generates 'k' sources from the 'belief' prior. """
        prior = np.copy(self.belief_) / self.belief_.sum()

        # Choose a 'k' random numbers uniformly in [0, 1]. These will be treated
        # as evaluations of the CDF, and since they are uniform in [0, 1], the
        # points at which they occur are distributed according to 'prior'.
        cdf_evals = sorted(np.random.uniform(0.0, 1.0, self.k_))

        # Walk the 'prior' until we get to each 'cdf_eval' and record which voxel
        # we are in.
        sources = []
        current_cdf_index = 0

        cdf = 0.0
        for ii in range(prior.shape[0]):
            for jj in range(prior.shape[1]):
                cdf += prior[ii, jj]

                # Check if we just passed the next 'cdf_eval'.
                if cdf >= cdf_evals[current_cdf_index]:
                    # Generate a new source here.
                    source = Source2D(float(ii) + 0.5, float(jj) + 0.5)
                    sources.append(source)

                    # Return if we've got enough sources.
                    if len(sources) == self.k_:
                        return sources

                    # Increment 'current_cdf_index'.
                    current_cdf_index += 1

        # Should never get here.
        assert(False)
        return []

    def Update(self, sensor, solve=True):
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

        # Identify voxels that are in range. Store row-major flattened indices.
        in_view_list = []
        for ii in range(self.belief_.shape[0]):
            for jj in range(self.belief_.shape[1]):
                if sensor.VoxelInView(ii, jj):
                    in_view_list.append(ii*self.belief_.shape[1] + jj)

        # Set up constrained least squares problem and solve.
        self.viewed_lists_.append(in_view_list)
        self.viewed_sets_.append(ImmutableSet(in_view_list))
        self.measurements_.append(measurement)

        if solve:
            self.SolveLeastSquares()

        return True

    def SolveLeastSquares(self):
        """
        Solve the following constrained linear least squares problem:
                 min 0.5 * sum_i ((sum_{j in S_i} p_j) - z_i)^2
                     + 0.5 * alpha * ((sum_j p_j) - k)^2
                 s.t. all probabilities are between 0 and 1
        where S_i = {j : voxel v_j is in view at time i}.

        This is a linear objective, but the constraints make it slightly
        annoying to solve in closed form, so for now we solve with non-linear
        least squares optimization, providing derivatives wherever possible.
        """
        belief = self.belief_.flatten()

        # Pre-compute the (constant) Jacobian matrix.
        J = np.empty((len(self.viewed_sets_) + 1, len(belief)), dtype=np.float)
        J[:-1, :] = map(lambda ii : map(lambda jj : (jj in self.viewed_sets_[ii]),
                                        range(len(belief))),
                        range(len(self.viewed_sets_)))
        J[-1:, :] = math.sqrt(self.alpha_)

        # Compute the vector of residuals.
        def residuals(p, viewed_lists, measurements, k, alpha, jac):
            resids = [(p[S].sum() - z) for (S, z) in zip(viewed_lists, measurements)]
            resids.append(math.sqrt(alpha) * (p.sum() - k))
            return np.array(resids)

        # Compute the Jacobian matrix.
        def jacobian(p, viewed_lists, measurements, k, alpha, jac):
            return jac

        try:
            result = least_squares(residuals, belief, jac=jacobian, bounds=(0.0, 1.0),
                                   args=(self.viewed_lists_, self.measurements_,
                                         self.k_, self.alpha_, J),
                                   xtol=0.001,
                                   verbose=0)
            self.belief_ = np.reshape(result.x, self.belief_.shape)

        except Exception, e:
            pass


    def SimulateTrajectory(self, sensor, trajectory, niters=1):
        """
        Return expected map entropy after taking scans at each GridPose2D
        in the trajectory.

        Expectation is based on Monte Carlo simulation using the specified
        number of iterations.
        """
        entropy_total = 0

        # Save the current state. This could be speeded up.
        current_belief = np.copy(self.belief_)
        current_viewed_lists = list(self.viewed_lists_)
        current_viewed_sets = list(self.viewed_sets_)
        current_measurements = list(self.measurements_)

        # Monte Carlo simulation.
        for ii in range(niters):
            for pose in trajectory:
                # Set sensor pose.
                sensor.ResetPose(pose)

                # Simulate sensor measurement.
                assert self.Update(sensor, False)

            # Once all poses have been simulated, compute entropy.
            self.SolveLeastSquares()
            entropy_total += self.Entropy()

            # Restore state.
            self.belief_ = np.copy(current_belief)
            self.viewed_lists_ = list(current_viewed_lists)
            self.viewed_sets_ = list(current_viewed_sets)
            self.measurements_ = list(current_measurements)

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
