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

#ifndef __CYCLOSTATIONARY_TYPES_H__
#define __CYCLOSTATIONARY_TYPES_H__

#include <cmath>
#include <complex>
#include <map>
#include <vector>

namespace cyclostationary {

#define MYPRINT false
#define myprintf(...)    \
  if (MYPRINT) {         \
    printf(__VA_ARGS__); \
  }
#define Vec(a, b) std::vector<__typeof(*(a))>((a), (a) + (b))

#define USE_FLOAT

#ifdef USE_FLOAT
typedef float real_t;
typedef std::complex<real_t> complex_t;
#endif

#ifdef USE_DOUBLE
typedef double real_t;
typedef std::complex<real_t> complex_t;
#endif

}  // namespace cyclostationary
#endif  // __CYCLOSTATIONARY_TYPES_H__
