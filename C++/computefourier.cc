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

#include "computefourier.h"


#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unordered_map>

#include "filters.h"
#include "types.h"
#include "utils.h"
namespace cyclostationary {

const int TimeSMod(const int& x, const int& a, const int& n) {
  return int((((long long int)x) * a) % n);
}

std::unique_ptr<std::vector<int>> Comb_Filt21(
    const std::vector<complex_t>& input_signal, int size_of_n, const int num,
    const int size_of_comb_window, std::vector<int>& indexn2) {
  if (size_of_n % size_of_comb_window) {
    fprintf(stderr, "Warning: W_Comb is not divisible by N, which algorithm expects.\n");
    assert(size_of_n % size_of_comb_window == 0);
  }

  std::vector<complex_t> vec_x_temp_samples(size_of_comb_window);

  // indexn2 store the relavent index from input_signal to x
  // (sigma * idx + offset) % M; sigma = size_of_n / size_of_comb_window;
  for (int i = 0; i < size_of_comb_window; i++) {
    vec_x_temp_samples[i] = input_signal[indexn2[i]];
  }

  // Compute M-point FFT; M = size_of_comb_window;
  fftw_dft(vec_x_temp_samples, size_of_comb_window, vec_x_temp_samples);

  std::vector<real_t> vec_samples(size_of_comb_window);
  for (int i = 0; i < size_of_comb_window; i++)
    vec_samples[i] = cabs2(vec_x_temp_samples[i]);

  // comb_approved: from M largest elements index of FFT(y) find largest 2k
  // elements.
  return find_largest_indices(num, vec_samples);
}

std::unique_ptr<std::vector<int>> inner_loop_locate1(
    const std::vector<complex_t>& input_signal, int total_size_n,
    std::vector<int>& indexn, const Filter& filter, const int num,
    const int size_of_bucket, std::vector<complex_t>& output_x_samp) {
  if (total_size_n % size_of_bucket) {
    fprintf(stderr, "Warning: W_Comb is not divisible by N, which algorithm expects.\n");
    assert(total_size_n % size_of_bucket == 0);
  }

  // size of x_sampt is B; B-point FFT; initial to 0;
  // complex_t* x_sampt = (complex_t*)malloc(B * sizeof(*x_sampt));
  std::vector<complex_t> x_samp_temp(size_of_bucket);

  // Permute, dot, collate all in one loop.
  // indexn store the matched index for filter multiplication
  // indexn[i] = (i * sigma + offset) % total_size_n; sigma = a_inverse;
  for (int i = 0; i < (int)filter.time.size(); i++) {
    // in inner loop, x_sampt = y = G*(P_{sigma,tau}x)
    x_samp_temp[i % size_of_bucket] += input_signal[indexn[i]] * filter.time[i];
  }

  // Step3:z^hat = y^hat_i(n/size_of_bucket)
  fftw_dft(output_x_samp, size_of_bucket, x_samp_temp);

  std::vector<real_t> norm_x_samp_temp(size_of_bucket);
  for (int i = 0; i < size_of_bucket; i++)
    norm_x_samp_temp[i] = cabs2(output_x_samp[i]);

  return find_largest_indices(num, norm_x_samp_temp);
}

int inner_loop_filter_regular(std::vector<int>& J_large_indices,
                              const int total_size_n, const int num,
                              const int size_of_bucket, const int sigma,
                              const int loop_threshold, std::vector<int>& score,
                              std::vector<int>& hits) {
  // Given the set of large samples, find the locations in [n] that map there
  // and output them
  int n = total_size_n;
  for (int i = 0; i < num; i++) {
    int low, high;
    low = (int(ceil((J_large_indices[i] - 0.5) * n / size_of_bucket)) + n) % n;
    high = (int(ceil((J_large_indices[i] + 0.5) * n / size_of_bucket)) + n) % n;
    int loc = TimeSMod(low, sigma, n);
    for (int j = low; j != high; j = (j + 1) % n) {
      score[loc]++;
      if (score[loc] == loop_threshold) hits.push_back(loc);
      loc = (loc + sigma) % n;
    }
  }

  return 0;
}

int inner_loop_filter_Comb(const std::vector<int>& J_large_indices,
                           const int total_size_n, const int num,
                           const int size_of_bucket, const int sigma,
                           const int sigma_inverse, const int loop_threshold,
                           std::vector<int>& score, std::vector<int>& hits,
                           const std::vector<int>& comb_approved,
                           const int size_of_comb_window) {
  int n = total_size_n;
  int num_Comb = comb_approved.size();
  // std::pair<int, int>* permuted_approved =
  // (__typeof(permuted_approved))malloc(
  //     num_Comb * sizeof(*permuted_approved));
  std::vector<std::pair<int, int>> permuted_approved;

  for (int m = 0; m < num_Comb; m++) {
    //(comb_approved[m]*ai) % size_of_comb_window
    int prev = TimeSMod(comb_approved[m], sigma_inverse, size_of_comb_window);
    // pair (comb_approved[m]*ai) % size_of_comb_window with (prev*a) % n
    permuted_approved.push_back(std::make_pair(prev, TimeSMod(prev, sigma, n)));
  }
  // sort it from small to large based on comb_index
  std::sort(permuted_approved.begin(), permuted_approved.end());

  // compute intersection of permuted_approved and indices close to
  // J_large_indices * n / B, then invert to get true locations.
  for (int i = 0; i < num; i++) {
    int low, high;
    low = (int(ceil((J_large_indices[i] - 0.5) * n / size_of_bucket)) + n) % n;
    high = (int(ceil((J_large_indices[i] + 0.5) * n / size_of_bucket)) + n) % n;
    // find the first permuted_approved.first(prev) larger than
    // (low%size_of_comb_window, -1) and compute the relative index of
    // permuted_approved
    int index =
        int(std::upper_bound(permuted_approved.begin(), permuted_approved.end(),
                             std::make_pair(low % size_of_comb_window, -1)) -
            permuted_approved.begin());
    // difference of index % M = 0
    int location = low - (low % size_of_comb_window);
    // the inverse location info: (location*a) % n
    int locinv = TimeSMod(location, sigma, n);
    // find the related location in permuted_approved
    for (int j = index;; j++) {
      // searching
      if (j == num_Comb) {
        // back to start position
        j = 0;
        // (location+M)%n, search every M point; M is size_of_comb_window
        location = (location + size_of_comb_window) % n;
        // update the inverse location
        locinv = TimeSMod(location, sigma, n);
      }
      // first: (comb_approved[m]*ai)%size_of_comb_window
      int approved_loc = location + permuted_approved[j].first;
      if ((low < high && (approved_loc >= high || approved_loc < low)) ||
          (low > high && (approved_loc >= high && approved_loc < low)))
        break;
      // second: (((comb_approved[m]*ai)%size_of_comb_window)*a)%n
      int loc = (locinv + permuted_approved[j].second) % n;
      score[loc]++;
      // loops for the number of location_loops
      // the index shows up more than loos_threshold will count
      if (score[loc] == loop_threshold) hits.push_back(loc);
    }
  }

  return 0;
}

std::map<int, complex_t> estimate_values(
    const std::vector<int>& hits, std::vector<std::vector<complex_t>>& x_samp,
    const int& total_loops, const int size_of_n,
    const std::vector<int>& permute, const int B_location,
    const int B2_estimation, const Filter& filter_location,
    const Filter& filter_estimation, const int location_loops) {
  int hits_found = hits.size();
  int n = size_of_n;
  std::map<int, complex_t> ans_output;

  for (int i = 0; i < hits_found; i++) {
    // save the candidate coeff and then estimate the target
    std::vector<std::vector<real_t>> values(2);

    for (int j = 0; j < total_loops; j++) {
      int cur_B = (j < location_loops) ? B_location : B2_estimation;
      const Filter* cur_filter =
          (j < location_loops) ? &filter_location : &filter_estimation;
      //(permute[j]*hits[i])%n
      int permuted_index = TimeSMod(permute[j], hits[i], n);
      int hashed_to = permuted_index / (n / cur_B);
      int dist = permuted_index % (n / cur_B);
      if (dist > (n / cur_B) / 2) {
        hashed_to = (hashed_to + 1) % cur_B;
        dist -= n / cur_B;
      }
      dist = (n - dist) % n;
      complex_t filter_value = cur_filter->freq[dist];
      complex_t temp = x_samp[j][hashed_to] / filter_value;
      values[0].push_back(temp.real());
      values[1].push_back(temp.imag());
    }

    int location = (total_loops - 1) / 2;
    // position save the length of values real and imag
    if (values[0].size() != values[1].size()) {
      fprintf(stderr, "Warning: Length of real part and imag part is not matched.\n");
      assert(0);
    }
    // get the median value, only the left values of location is <= location
    // value
    for (int a = 0; a < 2; a++)
      std::nth_element(values[a].begin(), values[a].begin() + location,
                       values[a].end());
    real_t realv = values[0][location];
    real_t imagv = values[1][location];
    ans_output[hits[i]] = complex_t(realv, imagv);
  }

  return ans_output;
}

std::map<int, complex_t> outer_loop21(
    const std::vector<complex_t>& input_signal,
    const SfftFuncParams sfft_func_parameter) {
  const int n = sfft_func_parameter.size_of_n;
  int total_loops = sfft_func_parameter.parameters.loc_loops +
                    sfft_func_parameter.parameters.est_loops;
  const std::vector<int>& indexn2 = *(sfft_func_parameter.index_trans);
  int B2_estimation = sfft_func_parameter.b_est;
  int num = sfft_func_parameter.b_thresh;
  int B_location = sfft_func_parameter.b_loc;
  int size_of_comb_window = sfft_func_parameter.size_of_comb_window;
  int Comb_loops = sfft_func_parameter.parameters.comb_loops;
  int loop_threshold = sfft_func_parameter.parameters.threshold_loops;
  int location_loops = sfft_func_parameter.parameters.loc_loops;
  const std::vector<int>& randnseed = *(sfft_func_parameter.sigma_randn);
  const Filter& filter_estimation = *(sfft_func_parameter.filter_estimation);
  const Filter& filter_location = *(sfft_func_parameter.filter_location);

  // for random pick the sigma for permutation
  std::vector<int> permute_sigma(total_loops);
  std::vector<int> permute_tau(total_loops);

  // estimate loops number (location + estimation) of vector
  std::vector<std::vector<complex_t>> x_samp(total_loops);

  for (int i = 0; i < total_loops; i++) {
    // IF: to save the result from location loops
    // ELSE: to save the result from estimation loops
    if (i < location_loops)
      x_samp[i].resize(B_location);
    else
      x_samp[i].resize(B2_estimation);
  }

  // score is large size of n
  std::vector<int> score(n);

  // hits in size of n
  std::vector<int> hits;

  // BEGIN Comb
  std::vector<int> comb_approved;
  int num_Comb = num;

  if (sfft_func_parameter.myversion.with_comb_filter) {
    for (int i = 0; i < Comb_loops; i++) {
      std::vector<int> tempindex;
      tempindex.insert(tempindex.begin(),
                       indexn2.begin() + i * size_of_comb_window,
                       indexn2.begin() + (i + 1) * size_of_comb_window);
      std::unique_ptr<std::vector<int>> temp_comb_approved =
          Comb_Filt21(input_signal, n, num, size_of_comb_window, tempindex);
      // comb_approved get target num largest coefficient of
      std::vector<int>& vec_temp = *temp_comb_approved.get();
      comb_approved.insert(comb_approved.end(), vec_temp.begin(),
                           vec_temp.end());
    }
  }

  if (Comb_loops > 1) {
    // sort comb_approved; comb_loops*num is size of comb_approved
    radix_sort(comb_approved, Comb_loops * num);
    int Last = 0;
    for (int i = 1; i < Comb_loops * num; i++) {
      // save unique values
      if (comb_approved[i] != comb_approved[Last])
        comb_approved[++Last] = comb_approved[i];
    }
    num_Comb = Last;
    comb_approved.erase(comb_approved.begin() + num_Comb, comb_approved.end());
  }

  // BEGIN INNER LOOPS
  for (int i = 0; i < total_loops; i++) {
    int sigma = randnseed[i + Comb_loops];
    // random() % n; 0 is fine; todo: set when modify autossca.cpp
    int offset = 0;
    // a*ai = 1 mod n
    int sigma_inverse = mod_inverse(sigma, n);
    permute_sigma[i] = sigma_inverse;
    permute_tau[i] = offset;

    int perform_location = (i < location_loops);
    // location loop need algorithm1, estimation loop no need for it

    // match the filter and bucket number for location loop or estimation loop
    const Filter* cur_filter =
        perform_location ? &filter_location : &filter_estimation;

    int cur_B = perform_location ? B_location : B2_estimation;
    std::vector<int> indexn;
    if (perform_location) {
      indexn.insert(indexn.begin(),
                    indexn2.begin() + Comb_loops * size_of_comb_window +
                        i * filter_location.time.size(),
                    indexn2.begin() + Comb_loops * size_of_comb_window +
                        (i + 1) * filter_location.time.size());
    } else {
      indexn.insert(
          indexn.begin(),
          indexn2.begin() + Comb_loops * size_of_comb_window +
              location_loops * filter_location.time.size() +
              (i - location_loops) * filter_estimation.time.size(),
          indexn2.begin() + Comb_loops * size_of_comb_window +
              location_loops * filter_location.time.size() +
              (i - location_loops + 1) * filter_estimation.time.size());
    }

    // x_samp save coeff values,J save num largest indices
    std::unique_ptr<std::vector<int>> J_large_coeff = inner_loop_locate1(
        input_signal, n, indexn, *cur_filter, num, cur_B, x_samp[i]);
    std::vector<int>& J_vec = *J_large_coeff.get();
    if (perform_location) {
      if (!sfft_func_parameter.myversion.with_comb_filter) {
        inner_loop_filter_regular(J_vec, n, num, cur_B, sigma, loop_threshold,
                                  score, hits);
      } else {
        // output score,hits
        inner_loop_filter_Comb(J_vec, n, num, cur_B, sigma, sigma_inverse,
                               loop_threshold, score, hits, comb_approved,
                               size_of_comb_window);
      }
    }
  }

  std::map<int, complex_t> ans = estimate_values(
      hits, x_samp, total_loops, n, permute_sigma, B_location, B2_estimation,
      filter_location, filter_estimation, location_loops);

  return ans;
}

}  // namespace cyclostationary