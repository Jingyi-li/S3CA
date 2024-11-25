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

#ifndef __CYCLOSTATIONARY_COMPUTEFOURIER_H__
#define __CYCLOSTATIONARY_COMPUTEFOURIER_H__

#include <complex>
#include <map>
#include <unordered_map>

#include "fftw.h"
#include "filters.h"
#include "parameters.h"
#include "types.h"

namespace cyclostationary {
typedef struct {
  bool with_comb_filter;
  int seed;
} SfftVersion;

typedef struct {
  SfftVersion myversion;
  const std::vector<int> *index_trans;
  int size_of_n;
  int np_channels;
  const std::vector<int> *sigma_randn;
  int kappa;
  const Filter *filter_location;
  const Filter *filter_estimation;
  int b_est;
  int b_thresh;
  int b_loc;
  int size_of_comb_window;
  SfftParameters parameters;
} SfftFuncParams;

const int TimeSMod(const int &x, const int &a, const int &n);

std::unique_ptr<std::vector<int>> Comb_Filt21(
    const std::vector<complex_t> &input_signal, int size_of_n, const int num,
    const int size_of_comb_window, std::vector<int> &indexn2);

/*
  Inner loop of the algorithm, part one.

  n-dimensional origx
  permute the fourier spectrum, take the first w coordinates.
  dot with the filter
  B-dimensional FFT
  return the top num samples.

  Output to
 */
std::unique_ptr<std::vector<int>> inner_loop_locate1(
    const std::vector<complex_t> &input_signal, int total_size_n,
    std::vector<int> &indexn, const Filter &filter, const int num,
    const int size_of_bucket, std::vector<complex_t> &output_x_samp);

/*
 Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after
 permutation.

 For each such i, increment score[i] and append to hits if score[i] reaches
 loop_threshold.
*/
int inner_loop_filter_regular(std::vector<int> &J_large_indices,
                              const int total_size_n, const int num,
                              const int size_of_bucket, const int sigma,
                              const int loop_threshold, std::vector<int> &score,
                              std::vector<int> &hits);

/*
 Find indices that (1) map to J_large_indices under the permutation and (2) lie
 in comb_approved mod size_of_comb_window.

 For each such i, increment hits[i] and append to hits_found if hits[i] reaches
 loop_threshold.

  */
int inner_loop_filter_Comb(const std::vector<int> &J_large_indices,
                           const int total_size_n, const int num,
                           const int size_of_bucket, const int sigma,
                           const int sigma_inverse, const int loop_threshold,
                           std::vector<int> &score, std::vector<int> &hits,
                           const std::vector<int> &comb_approved,
                           const int size_of_comb_window);
/*
  hits contains the indices that we want to estimate.

  x_samp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of x_samp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of x_samp.
 */
std::map<int, complex_t> estimate_values(
    const std::vector<int> &hits, std::vector<std::vector<complex_t>> &x_samp,
    const int &total_loops, const int size_of_n,
    const std::vector<int> &permute, const int B_location,
    const int B2_estimation, const Filter &filter_location,
    const Filter &filter_estimation, const int location_loops);
/*
  Outer loop of the algorithm.

  If we are performing the Comb heuristic, first we do so.

  Then, `loops` times:
    choose a random permutation
    run inner_loop_locate
    if in the first location_loops loops, also run inner_loop_filter

  at the end, `hits` contains the coordinates that appear at least
  loop_threshold of location_loops times.  We estimate the values at
  these coordinates as the median of the images x_samp[loops].

  Returns a map from coordinates to estimates.
 */

std::map<int, complex_t> outer_loop21(
    const std::vector<complex_t> &input_signal,
    const SfftFuncParams sfft_func_parameter);

}  // namespace cyclostationary
#endif  // __CYCLOSTATIONARY_COMPUTEFOURIER_H__
