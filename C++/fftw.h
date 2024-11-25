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
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#ifndef __CYCLOSTATIONARY_FFTW_H__
#define __CYCLOSTATIONARY_FFTW_H__

#include <vector>

#include "types.h"
#include<fftw3.h>

namespace cyclostationary {

#ifdef USE_FLOAT

#define fftw_plan_dft_1d fftwf_plan_dft_1d
#define fftw_plan fftwf_plan
#define fftw_execute fftwf_execute
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_free fftwf_free

#endif

// size_of_fft point FFT; backwards = 0 do forward fft; = 1 do backward fft;
int fftw_dft(std::vector<complex_t>& vec_output, int size_of_fft,
             std::vector<complex_t>& vec_input, int backwards = 0);

}  // namespace cyclostationary
#endif  // __CYCLOSTATIONARY_FFTW_H__