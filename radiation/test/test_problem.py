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
# Unit tests for Sensor2D class.
#
###########################################################################

import numpy as np
import numpy.linalg as LA
import math
from numpy.testing import assert_equal, assert_array_less

from sensor_2d import Sensor2D
from source_2d import Source2D
from grid_pose_2d import GridPose2D
from problem import Problem

""" Test that estimated distributions converge with more samples. """
def test_distribution_convergence():
    # Max deviation at any point in estimated matrices/vectors.
    kPrecision = 0.05

    # Define hyperparameters.
    kNumSamples = 10000
    kNumRows = 5
    kNumCols = 5
    kNumSources = 1
    kNumSteps = 2
    kAngularStep = 0.25 * math.pi
    kSensorParams = {"x" : kNumRows/2,
                     "y" : kNumCols/2,
                     "angle" : 0.0,
                     "fov" : 0.5 * math.pi}

    # Generate conditionals from the specified pose.
    pose = GridPose2D(kNumRows, kNumCols,
                      int(np.random.uniform(0.0, kNumRows)) + 0.5,
                      int(np.random.uniform(0.0, kNumCols)) + 0.5,
                      np.random.uniform(0.0, 2.0 * math.pi))

    # Create a problem.
    problem = Problem(kNumRows, kNumCols, kNumSources, kNumSteps,
                      kAngularStep, kSensorParams, kNumSamples)

    (pzx1, hzm1, traj_ids1) = problem.GenerateConditionals(pose)
    (pzx2, hzm2, traj_ids2) = problem.GenerateConditionals(pose)

    # Check that distributions are close.
    assert_equal(traj_ids1, traj_ids2)
    assert_array_less(abs(pzx1 - pzx2).max(), kPrecision)
    assert_array_less(abs(hzm1 - hzm2).max(), kPrecision)

    # Check that cost vectors are not too different.
    c1 = pzx1.T * hzm1
    c2 = pzx2.T * hzm2
    assert_array_less(abs(c1 - c2).max(), kPrecision)
