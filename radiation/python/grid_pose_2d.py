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
# Defines a pose (position + orientation) on a 2D grid of the specified
# dimensions.
#
###########################################################################

import numpy as np
import math

class GridPose2D:
    def __init__(self, nrows, ncols, x, y, angle):
        """ Constructor. """
        self.nrows_ = nrows
        self.ncols_ = ncols
        self.x_ = math.floor(x) + 0.5
        self.y_ = math.floor(y) + 0.5
        self.angle_ = angle

    @classmethod
    def Copy(self, pose):
        """ Initialize from an existing pose. """
        return GridPose2D(pose.nrows_, pose.ncols_, pose.x_, pose.y_, pose.angle_)

    def MoveBy(self, delta_x, delta_y, delta_angle):
        """ Move this pose by the specified delta in each dimension. """
        new_x = self.x_ + delta_x
        new_y = self.y_ + delta_y
        new_angle = self.angle_ + delta_angle

        # Catch going out of bounds.
        if not ((new_x >= 0) and (new_x < self.nrows_) and
                (new_y >= 0) and (new_y < self.ncols_)):
            return False

        # Not going out of bounds, so update coordinates.
        self.x_ = new_x
        self.y_ = new_y
        self.angle_ = new_angle
        return True
