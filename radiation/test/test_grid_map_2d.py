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
# Unit tests for GridMap2D class.
#
###########################################################################

import numpy as np
from numpy import linalg as LA
import math

from radiation.sensor_2d import Sensor2D
from radiation.source_2d import Source2D
from radiation.grid_map_2d import GridMap2D

""" Test that we can sense an empty map properly with a full FOV sensor. """
def test_empty_map():
    # Create a grid map, specifying a non-zero number of sources.
    kNrows = 10
    kNcols = 10
    kNsources = 1
    grid = GridMap2D(kNrows, kNcols, kNsources)

    # Initialize a sensor in the middle of the grid, with a 2pi FOV.
    kFieldOfView = 2.0 * math.pi
    params = {"x" : 0.5 * kNrows,
              "y" : 0.5 * kNcols,
              "fov" : kFieldOfView,
              "angle" : 0.0}
    sensor = Sensor2D(params, [])

    # Update a bunch of times, and check that grid has converged to zero.
    kNumUpdates = 10
    kEpsilon = 1e-4
    for ii in range(kNumUpdates):
        assert grid.Update(sensor)

    assert LA.norm(grid.belief_) < kEpsilon
