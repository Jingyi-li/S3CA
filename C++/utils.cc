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

#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>

#include "types.h"

namespace cyclostationary {

// 2^8
constexpr int kRadixCount = 256;
// mask a byte
constexpr int kByteMask = 0xFF;
// this is in byte and each byte is 8 bits; radix is sort in bits
constexpr int kRadixByte2Bits = 3;

// x[i] <- x[i-r]
void shift(std::vector<complex_t>& x, int r) {
  int n = x.size();
  r = (n + r) % n;
  assert(n >= r);
  std::unique_ptr<std::vector<complex_t>> temp_vector =
      std::make_unique<std::vector<complex_t>>(x.begin() + n - r, x.end());
  std::memmove(x.data() + r, x.data(), (n - r) * sizeof(complex_t));
  std::memcpy(x.data(), temp_vector->data(), r * sizeof(complex_t));
}

// Compute the greatest common divisor of a and b
// assumes a, b > 0.
int gcd(int a, int b) {
  if (a % b == 0) return b;
  return gcd(b, a % b);
}

real_t phase(complex_t x) { return std::arg(x); }

real_t cabs2(complex_t x) { return std::norm(x); }

// crappy inversion code I stole from elsewhere
// Undefined if gcd(a, n) > 1
int mod_inverse(int a, int n) {
  int i = n, v = 0, d = 1;
  while (a > 0) {
    // save the result in t and residual in a;
    int t = i / a, x = a;
    a = i % x;
    i = x;
    x = d;
    d = v - t * x;
    v = x;
  }
  v %= n;
  if (v < 0) v = (v + n) % n;
  return v;
}

/*
  Compute the num'th smallest element of the length n input.
  uses std::nth_element, but doesn't mutate input.
*/
real_t nth_element_immutable(const std::vector<real_t>& input, int num) {
  std::unique_ptr<std::vector<real_t>> temp_vector =
      std::make_unique<std::vector<real_t>>(input.begin(), input.end());
  std::nth_element(temp_vector->begin(), temp_vector->begin() + num,
                   temp_vector->end());
  std::vector<real_t>& temp_vector_addr = *temp_vector.get();
  real_t ans = temp_vector_addr[num];
  return ans;
}

/*
  Output the indices corresponding to the num largest elements of samples.
  Output is sorted.
*/
std::unique_ptr<std::vector<int>> find_largest_indices(
    int num, const std::vector<real_t>& samples) {
  int n = samples.size();
  assert(n >= num + 1);
  std::unique_ptr<std::vector<int>> output =
      std::make_unique<std::vector<int>>();
  std::vector<int>& vec_output = *output.get();
  // use num+1 so we can use > cutoff and probably get exactly num.
  // if we get fewer, the second pass uses == cutoff.
  real_t cutoff = nth_element_immutable(samples, n - num - 1);

  // output save the location index
  for (int i = 0; i < n; i++)
    if (samples[i] > cutoff) {
      vec_output.push_back(i);
    }
  if (int(vec_output.size()) < num) {
    for (int i = 0; i < n; i++) {
      if (samples[i] == cutoff) {
        vec_output.push_back(i);
        if (int(vec_output.size()) >= num) break;
      }
    }
    std::sort(vec_output.begin(), vec_output.end());
  }
  assert(int(vec_output.size()) == num);
  return output;
}

void radix(int byte, int size, std::vector<int>& vec_input,
           std::vector<int>& vec_temp) {
  std::vector<int> count(kRadixCount);

  byte = byte << kRadixByte2Bits;

  for (int i = 0; i < size; ++i)
    ++count[((vec_input[i]) >> (byte)) & kByteMask];
  for (int i = 1; i < kRadixCount; ++i) count[i] += count[i - 1];
  for (int i = size - 1; i >= 0; --i) {
    vec_temp[count[(vec_input[i] >> (byte)) & kByteMask] - 1] = vec_input[i];
    --count[(vec_input[i] >> (byte)) & kByteMask];
  }
}

void radix_sort(std::vector<int>& vec_input, int size) {
  // save intermediate values
  std::unique_ptr<std::vector<int>> uptr_temp =
      std::make_unique<std::vector<int>>(size);
  std::vector<int>& vec_temp = *uptr_temp.get();

  // based on bits of each byte to sort
  for (unsigned int i = 0; i < sizeof(int); i += 2) {
    // even byte
    radix(i, size, vec_input, vec_temp);

    // odd byte
    radix(i + 1, size, vec_temp, vec_input);
  }
}

void radix_filt(int byte, int size, std::vector<int>& vec_input,
                std::vector<int>& vec_temp, std::vector<complex_t>& vec_filter,
                std::vector<complex_t>& vec_temp_filter) {
  std::vector<int> count(kRadixCount);

  byte = byte << kRadixByte2Bits;

  for (int i = 0; i < size; ++i)
    ++count[((vec_input[i]) >> (byte)) & kByteMask];
  for (int i = 1; i < kRadixCount; ++i) count[i] += count[i - 1];
  for (int i = size - 1; i >= 0; --i) {
    vec_temp[count[(vec_input[i] >> (byte)) & kByteMask] - 1] = vec_input[i];
    vec_temp_filter[count[(vec_input[i] >> (byte)) & kByteMask] - 1] =
        vec_filter[i];
    --count[(vec_input[i] >> (byte)) & kByteMask];
  }
}

void radix_sort_filt(std::vector<int>& vec_input,
                     std::vector<complex_t>& vec_filter, int size) {
  std::unique_ptr<std::vector<int>> uptr_temp =
      std::make_unique<std::vector<int>>(size);
  std::vector<int>& vec_temp = *uptr_temp.get();
  std::unique_ptr<std::vector<complex_t>> uptr_temp_filer =
      std::make_unique<std::vector<complex_t>>(size);
  std::vector<complex_t>& vec_temp_filter = *uptr_temp_filer.get();

  for (unsigned int i = 0; i < sizeof(int); i += 2) {
    // even byte
    radix_filt(i, size, vec_input, vec_temp, vec_filter, vec_temp_filter);

    // odd byte
    radix_filt(i + 1, size, vec_temp, vec_input, vec_temp_filter, vec_filter);
  }
}

int floor_to_pow2(double x) {
  unsigned int ans;
  for (ans = 1; ans <= x; ans <<= 1)
    ;
  return ans / 2;
}

double binomial_cdf(double prob, int n, int needed) {
  double ans = 0;
  double choose = 1;
  for (int i = n; i >= needed; i--) {
    ans += choose * pow(prob, i) * pow(1 - prob, n - i);
    choose = choose * i / (n - i + 1);
  }
  return ans;
}

}  // namespace cyclostationary