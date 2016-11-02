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
# Generate an instance of the LP formulation and save the data.
#
###########################################################################

from problem import Problem
from grid_pose_2d import GridPose2D

import numpy as np
import math

# Files to save to.
pzx_file = "pzx_5x5_1000.npy"
hmz_file = "hmz_5x5_1000.npy"

# Define hyperparameters.
kNumSamples = 1000
kNumRows = 5
kNumCols = 5
kNumSources = 1
kNumSteps = 1
kAngularStep = 0.25 * math.pi
kSensorParams = {"x" : kNumRows/2,
                 "y" : kNumCols/2,
                 "angle" : 0.0,
                 "fov" : 0.5 * math.pi}

# Create a problem.
problem = Problem(kNumRows, kNumCols, kNumSources, kNumSteps,
                  kAngularStep, kSensorParams, kNumSamples)

# Generate conditionals from the specified pose.
pose = GridPose2D(kNumRows, kNumCols, kNumRows/2, kNumCols/2, 0.0)
(pzx, hmz) = problem.GenerateConditionals(pose)

print "P_{Z|X} shape: " + str(pzx.shape)
print "h_{M|Z} shape: " + str(hmz.shape)

np.save(pzx_file, pzx)
np.save(hmz_file, hmz)
print "Successfully saved to disk."
