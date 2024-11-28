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

#include "filters.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "fftw.h"
#include "utils.h"

namespace cyclostationary {

// The Chebyshev polynomials for Dolph Chebyshev Window
// mth order
double Cheb(const double order_m, const double x) {
  if (std::fabs(x) <= 1)
    return std::cos(order_m * std::acos(x));
  else
    return std::real(std::cosh((std::complex<double>)order_m *
                               std::acosh((std::complex<double>)x)));
}

// size_of_w = int((1 / M_PI) * (1 / lobefrac) * acosh(1. / tolerance));
std::vector<complex_t> make_dolphchebyshev_t(double lobefrac, double tolerance,
                                             const int size_of_w) {
  // make sure that w is odd
  assert(size_of_w % 2);

  std::vector<complex_t> vec_output(size_of_w);

  double t0 = std::cosh(std::acosh(1 / tolerance) / (size_of_w - 1));
  for (int i = 0; i < size_of_w; i++) {
    vec_output[i] =
        Cheb(size_of_w - 1, t0 * std::cos(M_PI * i / size_of_w)) * tolerance;
  }
  fftw_dft(vec_output, size_of_w, vec_output);
  shift(vec_output, size_of_w / 2);
  for (int i = 0; i < size_of_w; i++) vec_output[i] = std::real(vec_output[i]);
  return vec_output;
}

double sinc(double x) {
  if (x == 0) return 1;
  return std::sin(x) / (x);
}

// compute the frequency domain and time domain of the filter based on given
// parameters
void make_multiple_t(Filter &filter, const int n, const int b) {
  int size_of_w = filter.time.size();
  assert(b <= n);
  assert(size_of_w <= n);
  std::vector<complex_t> vec_temp(n);
  filter.freq.resize(n);
  std::memcpy(vec_temp.data(), filter.time.data() + size_of_w / 2,
              (size_of_w - (size_of_w / 2)) * sizeof(complex_t));
  std::memcpy(vec_temp.data() + n - size_of_w / 2, filter.time.data(),
              (size_of_w / 2) * sizeof(complex_t));
  fftw_dft(vec_temp, n, vec_temp);
  complex_t tempsum = complex_t(0, 0);
  for (int i = 0; i < b; i++) tempsum += vec_temp[i];

  real_t max = 0;
  int offset = b / 2;
  for (int i = 0; i < n; i++) {
    filter.freq[(i + n + offset) % n] = tempsum;
    max = std::max(max, std::abs(tempsum));
    tempsum = tempsum + (vec_temp[(i + b) % n] - vec_temp[i]);
  }
  for (int i = 0; i < n; i++) filter.freq[i] /= max;

  complex_t offsetc = complex_t(1, 0);
  complex_t step = std::exp(complex_t(0, -2 * M_PI * (size_of_w / 2) / n));
  for (int i = 0; i < n; i++) {
    filter.freq[i] *= offsetc;
    offsetc *= step;
  }

  fftw_dft(vec_temp, n, filter.freq, 1);
  std::memmove(filter.time.data(), vec_temp.data(),
               size_of_w * sizeof(complex_t));

  for (int i = 0; i < size_of_w; i++) filter.time[i] /= n;
}
}  // namespace cyclostationary