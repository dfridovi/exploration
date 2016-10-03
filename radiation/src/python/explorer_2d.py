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

from grid_map_2d import GridMap2D
from source_2d import Source2D

import numpy as np

class Explorer2D:
    def __init__(self, nrows, ncols, k):
        """
        Constructor. Takes in dimensions and number of sources.
        Sets up a new GridMap2D to keep track of belief over time.
        Generates random sources.
        """
        self.map_ = GridMap2D(nrows, ncols, k)
        self.sources_ = []

        for ii in range(k):
            x = np.random.uniform(0.0, float(nrows))
            y = np.random.uniform(0.0, float(ncols))
            self.sources_.append(Source2D(x, y))

    def PlanAhead(self, nsteps, N, M):
        """
        Simulate N random nsteps trajectories, each using M iterations.
        Return the best one.
        """
        # TODO!
