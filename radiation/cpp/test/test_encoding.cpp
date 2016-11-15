/*
 * Copyright (c) 2016, The Regents of the University of California (Regents).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the name of the copyright holder nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Please contact the author(s) of this library if you have any questions.
 * Authors: David Fridovich-Keil   ( dfk@eecs.berkeley.edu )
 */

///////////////////////////////////////////////////////////////////////////////
//
// Unit tests for encoding/decoding sources, trajectories, and measurements.
//
///////////////////////////////////////////////////////////////////////////////

#include "../src/encoding.h"

#include <gtest/gtest.h>
#include <vector>
#include <random>

using namespace radiation;

// Test encoding/decoding for sources.
TEST(Encoding, TestSources) {
  const unsigned int kNumRows = 10;
  const unsigned int kNumCols = 10;
  const unsigend int kNumSources = 5;

  // Make a random number generator for each dimension.
  std::random_device rd;
  std::default_random_engine rng(rd);
  std::uniform_int_distribution<unsigned int> unif_rows(0, kNumRows - 1);
  std::uniform_int_distribution<unsigned int> unif_cols(0, kNumCols - 1);

  // Generate a buch of random sources.
  std::vector<Source2D> sources;
  for (size_t ii = 0; ii < kNumSources; ii++)
    sources.push_back(Source2D(unif_rows(rng), unif_cols(rng)));

  // Encode these sources.
  const unsigned int map_id = EncodeMap(sources, kNumRows, kNumCols);

  // Decode the id back into a list of sources.
  std::vector<Source2D> decoded_sources;
  DecodeMap(map_id, kNumRows, kNumCols, kNumSources, decoded_sources);
}
