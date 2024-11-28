/*
 * Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *   Eric Price, Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITfftw_dft
 *
 */

#include "fftw.h"

#include <array>
#include <map>
#include <vector>

namespace cyclostationary {

int fftw_dft(std::vector<complex_t>& vec_output, int size_of_fft,
             std::vector<complex_t>& vec_input, int backwards) {
  fftw_plan p;
  // Convert from complex<real_t> to 2D array of real_t
  real_t(*ptr_input)[2] = reinterpret_cast<real_t(*)[2]>(vec_input.data());
  real_t(*ptr_output)[2] = reinterpret_cast<real_t(*)[2]>(vec_output.data());
  p = fftw_plan_dft_1d(size_of_fft, ptr_input, ptr_output,
                       backwards ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);
  return 0;
}
}  // namespace cyclostationary