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
# Defines a 2D radiation sensor, which can sense gamma rays in a particular
# field of view (parameterized by direction and angular width), from
# a particular location.
#
###########################################################################

import numpy as np

class Sensor2D:
    def __init__(self, params):
        """
        Constructor. Parameters must include:
        1. the index of the grid cell on which the sensor is located
        2. the unit vector orientation of the sensor
        3. the angular field of view (in radians)
        """
        self.x_ = params["x"]
        self.y_ = params["y"]
        self.fov_ = params["fov"]
        self.orientation_ = params["orientation"]

    def Sense(self, sources):
        """
        Sense the given sources. Returns the number of sources in the
        field of view.
        """
        num_in_view = 0

        # Check each source, and see if it is in the field of view.
        for source in sources:
            if InView(source):
                num_in_view += 1

        return num_in_view

    def InView(self, source):
        """ Check if a source is in the field of view. """

        # Get vector to source.
        dx = source.x_ - self.x_
        dy = source.x_ - self.x_

        # In view if the angle of this vector in the plane is between
        # our orientation +/- half the field of view.
