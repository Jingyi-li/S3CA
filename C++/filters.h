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

#ifndef __CYCLOSTATIONARY_FILTERS_H__
#define __CYCLOSTATIONARY_FILTERS_H__

#include <memory>

#include "types.h"

namespace cyclostationary {

typedef struct {
  std::vector<complex_t> time;
  std::vector<complex_t> freq;
} Filter;

/*
  Create a window function such that:
      the main lobe has width 2 * n * filterfrac
      outside the main lobe, the linf residual is tolerance
  Computes the required w.
  Allocates and returns the filter.
 */
double Cheb(const double order_m, const double x);

// Generate dolphchebyshev window based on window_size lobefrac and tolerance
std::vector<complex_t> make_dolphchebyshev_t(double lobefrac, double tolerance,
                                             const int size_of_w);
double sinc(double x);

/*
  Modifies a w-dimensional window function to have n-dimensional FFT
  the sum of b adjacent ones previously.

  Allocates and returns a Filter instance pointing to the modified x and an
  n-dimensional FFT of it.
 */
void make_multiple_t(Filter &filter, const int n, const int b);

}  // namespace cyclostationary
#endif  // __CYCLOSTATIONARY_FILTERS_H__
