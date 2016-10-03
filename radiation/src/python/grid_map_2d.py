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

class GridMap2D:
    def __init__(self, nrows, ncols, k):
        """
        Constructor. Takes in dimensions and number of sources k, and
        and generates a uniform prior.
        """
        self.belief_ = np.ones((nrows, ncols), dtype=np.float) / k

    def Update(self, sensor_params, sensor_reading):
        """
        Update belief about the world, given that the specified sensor
        reading was received by a sensor with the specified parameters.
        """
        # TODO!

    def Simulate(self, sensor_params, niters):
        """
        Return expected map entropy after receiving a measurement from
        the specified location/orientation. Expectation is based on
        Monte Carlo simulation using the specified number of iterations.
        """
        # TODO!

    def Entropy(self):
        """ Compute the entropy of the map. """
        # TODO!
