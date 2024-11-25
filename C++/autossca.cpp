#include "autossca.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <map>
#include <random>
#include <unordered_map>

#include "SSCA.h"
#include "computefourier.h"
#include "fftw.h"
#include "filters.h"
#include "parameters.h"
#include "types.h"
#include "utils.h"

namespace cyclostationary {

Filter make_filter(int size_of_n, double Bcst, int kappa, double tolerance,
                   double factor) {
  // define the parameters of filter window
  real_t BB =
      (unsigned)(Bcst * sqrt((real_t)size_of_n * kappa / (log2(size_of_n))));
  real_t lobefrac = 0.5 / (BB);
  int b = int(factor * ((real_t)size_of_n / BB));

  Filter filter;
  int size_of_w = int((1 / M_PI) * (1 / lobefrac) * acosh(1. / tolerance));
  // make sure that w is odd
  if (!(size_of_w % 2)) size_of_w--;

  filter.time = make_dolphchebyshev_t(lobefrac, tolerance, size_of_w);
  make_multiple_t(filter, size_of_n, b);
  return filter;
}

// Compute the location information inside of SFFT with relative position
std::map<unsigned int, real_t> SFFTfunc33(
    std::vector<complex_t> &input_signal,
    const SfftFuncParams sfft_func_parameters, int indx_channel) {
  std::map<int, complex_t> ans;
  std::map<unsigned int, real_t> Sx;
  ans = outer_loop21(input_signal, sfft_func_parameters);

  int num_candidates = (int)ans.size();
  std::vector<std::pair<real_t, int>> candidates;
  int kappa = sfft_func_parameters.kappa,
      size_of_n = sfft_func_parameters.size_of_n,
      np_channels = sfft_func_parameters.np_channels;

  for (__typeof(ans.begin()) it = ans.begin(); it != ans.end(); it++) {
    int key = it->first;
    complex_t value = it->second;
    candidates.push_back(std::make_pair(std::abs(value), key));
  }
  /*for some case the number of num_candidates is less than k*/
  int k2 = (kappa >= num_candidates) ? (int)((num_candidates)*3 / 4) : kappa;
  std::nth_element(candidates.begin(), candidates.begin() + num_candidates - k2,
                   candidates.end());
  for (int i = 0; i < k2; i++) {
    int key = candidates[num_candidates - k2 + i].second;
    // We are more interested about the absolute values
    // return the absolute value of sFFT and map to SCD matrix
    key = (key + size_of_n / 2) % size_of_n;
    float alpha =
        ((float)key) / size_of_n + ((float)indx_channel) / np_channels;
    float freq =
        (((float)indx_channel) / np_channels - ((float)key) / size_of_n) / 2;
    int idxfreq = round(np_channels * (freq + 0.5));
    int idxalpha = round(size_of_n * (alpha));
    Sx[idxfreq + idxalpha * np_channels] =
        candidates[num_candidates - k2 + i].first;
  }
  return Sx;
}

std::unique_ptr<std::vector<real_t>> autossca_fftw(
    const std::vector<complex_t> &input_signal, int np_channels,
    int size_of_n) {
  int i, j, k, m, idxfreq, idxalpha;
  std::unique_ptr<std::vector<real_t>> Sx_ptr =
      std::make_unique<std::vector<real_t>>();
  std::vector<real_t> &Sx = *Sx_ptr.get();
  real_t alpha, frequency;

  if (np_channels < 2) {
    printf("Error: np_channels must be at least 2.\n");
    assert(0);
  }

  if (size_of_n < np_channels) {
    printf("Error: size_of_n must be at least np_channels.\n");
    assert(0);
  }

  /*
  *******************
  Windowing
  *******************
  */
  std::vector<complex_t> XW(np_channels * size_of_n);
  for (i = 0; i < size_of_n; i++) {
    for (j = 0; j < np_channels; j++) {
      XW[j + i * np_channels] = input_signal[i + j] * chebwin_64[j];
    }
  }

  /*
  *******************
  First FFT
  *******************
  */

  std::vector<complex_t> input_XFx1(np_channels);
  std::vector<complex_t> output_XFf1(np_channels);
  real_t(*ptr_input)[2] = reinterpret_cast<real_t(*)[2]>(input_XFx1.data());
  real_t(*ptr_output)[2] = reinterpret_cast<real_t(*)[2]>(output_XFf1.data());

  fftw_plan p;

  p = fftw_plan_dft_1d(np_channels, ptr_input, ptr_output, FFTW_FORWARD,
                       FFTW_MEASURE);
  for (j = 0; j < size_of_n; j++) {
    for (i = 0; i < np_channels; i++) {
      input_XFx1[i] = XW[i + j * np_channels];
    }
    fftw_execute(p);
    // with frequency shift
    for (i = 0; i < np_channels / 2; i++) {
      XW[i + j * np_channels] = output_XFf1[i + np_channels / 2];
    }
    // with frequency shift
    for (i = np_channels / 2; i < np_channels; i++) {
      XW[i + j * np_channels] = output_XFf1[i - np_channels / 2];
    }
  }
  fftw_destroy_plan(p);

  /*
*******************************************
Down Conversion + Conjugate Multiplication
*******************************************
*/

  std::vector<complex_t> Et;

  for (k = 0; k < np_channels; k++) {
    real_t coeff_idx = real_t(k) / np_channels;
    Et.push_back(
        complex_t(cos(-2.0 * M_PI * coeff_idx), sin(-2.0 * M_PI * coeff_idx)));
  }

  for (m = 0; m < size_of_n; m++) {
    complex_t xc = std::conj(input_signal[m + np_channels / 2]);
    for (k = 0; k < np_channels; k++) {
      XW[k + m * np_channels] =
          Et[(((k + np_channels / 2) % np_channels) * m) % np_channels] *
          XW[k + m * np_channels] * xc;
    }
  }

  /*
  ************************
  Second FFT
  ************************
  */

  std::vector<complex_t> input_XFx2(size_of_n);
  std::vector<complex_t> output_XFf2(size_of_n);
  real_t(*ptr_input2)[2] = reinterpret_cast<real_t(*)[2]>(input_XFx2.data());
  real_t(*ptr_output2)[2] = reinterpret_cast<real_t(*)[2]>(output_XFf2.data());

  fftw_plan p2;
  p2 = fftw_plan_dft_1d(size_of_n, ptr_input2, ptr_output2, FFTW_FORWARD,
                        FFTW_ESTIMATE);
  for (j = 0; j < np_channels; j++) {
    for (i = 0; i < size_of_n; i++) {
      input_XFx2[i] = XW[i * np_channels + j];
    }

    fftw_execute(p2);

    // with frequency shift
    for (i = 0; i < size_of_n / 2; i++) {
      XW[i * np_channels + j] =
          output_XFf2[i + size_of_n / 2] / complex_t(size_of_n, 0);
    }
    for (i = size_of_n / 2; i < size_of_n; i++) {
      XW[i * np_channels + j] =
          output_XFf2[i - size_of_n / 2] / complex_t(size_of_n, 0);
    }
  }
  fftw_destroy_plan(p2);

  /*
  ************************
  SCD Matrix
  ************************
  */
  double alpha_start = -1.0;
  double alpha_step = 1.0 / size_of_n;
  Sx.resize(2 * size_of_n * np_channels);

  std::vector<real_t> alphao(2 * size_of_n);
  std::vector<real_t> freqo(np_channels);
  for (i = 0; i < 2 * size_of_n; i++) alphao[i] = alpha_start + i * alpha_step;
  double freq_start = -0.5;
  double freq_step = 1.0 / np_channels;
  for (i = 0; i < np_channels; i++) freqo[i] = freq_start + i * freq_step;

  for (i = 0; i < size_of_n; i++) {
    for (j = 0; j < np_channels; j++) {
      alpha = ((float)i) / size_of_n + ((float)j) / np_channels;
      frequency = (((float)j) / np_channels - ((float)i) / size_of_n) / 2;
      idxfreq = round(np_channels * (frequency + 0.5));
      idxalpha = round(size_of_n * (alpha));
      Sx[idxfreq + idxalpha * np_channels] = std::abs(XW[i * np_channels + j]);
    }
  }
  return Sx_ptr;
}

std::vector<int> GenerateSigmaRandom(std::mt19937 &generator, int size_of_n,
                                     int size_of_comb, int comb_loops,
                                     int sumloops) {
  std::vector<int> sigma_randn;
  // for comb loops the offset = random_value % size_of_comb
  std::uniform_int_distribution<int> comb_distrib(0, size_of_comb - 1);
  // for location loops and estimation loops; the sigma should be co prime with
  // size_of_n in range of [0, size_of_n)]
  std::uniform_int_distribution<int> n_distrib(0, size_of_n - 1);

  for (int ii = 0; ii < sumloops; ii++) {
    if (ii < comb_loops) {
      // Comb filter offset this is random
      sigma_randn.push_back(comb_distrib(generator));

    } else {
      // until sigma is odd
      int new_value = n_distrib(generator);
      while (gcd(new_value, size_of_n) != 1) {
        new_value = n_distrib(generator);
      }
      sigma_randn.push_back(new_value);
    }
  }
  return sigma_randn;
}

std::map<unsigned int, real_t> sparsessca7(
    const SfftVersion &myversion, const std::vector<complex_t> &input_signal,
    int np_channels, int size_of_n, int kappa, const SfftParameters &parameters,
    const Filter &filter_location, const Filter &filter_estimation) {
  double Bcst_loc = parameters.loc_bcst, Bcst_est = parameters.est_bcst,
         Comb_cst = parameters.comb_cst;
  int comb_loops = parameters.comb_loops, loc_loops = parameters.loc_loops,
      est_loops = parameters.est_loops;

  int w_loc = filter_location.time.size();
  int w_est = filter_estimation.time.size();

  if (input_signal.size() == 0) {
    printf("Error: NULL input.\n");
    assert(0);
  }

  if (np_channels < 2) {
    printf("Error: np_channels must be at least 2.\n");
    assert(0);
  }

  if (size_of_n < np_channels) {
    printf("Error: size_of_n must be at least np_channels.\n");
    assert(0);
  }
  /*
  *******************
  Spare preparation
  *******************
  */

  size_of_n = floor_to_pow2(size_of_n);

  real_t BB_loc = (unsigned)(Bcst_loc * sqrt((real_t)size_of_n * kappa /
                                             (log2(size_of_n))));
  real_t BB_est = (unsigned)(Bcst_est * sqrt((real_t)size_of_n * kappa /
                                             (log2(size_of_n))));

  int B_loc = floor_to_pow2(BB_loc);
  int B_thresh = 2 * kappa;
  int B_est = floor_to_pow2(BB_est);

  int w_comb = 1;
  int sigma_comb = 1;
  if (myversion.with_comb_filter) {
    w_comb = floor_to_pow2(Comb_cst * size_of_n / B_loc);
    sigma_comb = size_of_n / w_comb;
  }

  /*Base on the loops and random seed to decide the seeds*/

  int sumloops = comb_loops + est_loops + loc_loops;
  int indexnsize = comb_loops * w_comb + w_loc * loc_loops + w_est * est_loops;
  std::vector<int> indexn_input(indexnsize), indexn_transform(indexnsize);

  std::mt19937 generator(myversion.seed * 171717);
  std::vector<int> sigma_randn = GenerateSigmaRandom(
      generator, size_of_n, sigma_comb, comb_loops, sumloops);

  for (int ii = 0; ii < sumloops; ii++) {
    if (ii < comb_loops) {
      for (int i = 0; i < w_comb; i++)
        indexn_input[i + ii * w_comb] = sigma_randn[ii] + i * sigma_comb;
    } else {
      // This should be offset but offset can be zero in the location and
      // estimation
      int index = 0;

      int per_loc = ii < comb_loops + loc_loops;
      int winloops = per_loc ? w_loc : w_est;
      int ai = mod_inverse(sigma_randn[ii], size_of_n);
      int startindex = per_loc
                           ? comb_loops * w_comb + (ii - comb_loops) * winloops
                           : comb_loops * w_comb + loc_loops * w_loc +
                                 (ii - comb_loops - loc_loops) * winloops;
      for (int i = 0; i < winloops; i++) {
        indexn_input[i + startindex] = index;
        index = (index + ai) % size_of_n;
      }
    }
  }

  // make the indexn_input unique
  std::memcpy(indexn_transform.data(), indexn_input.data(),
              indexnsize * sizeof(int));
  std::unordered_map<int, int> indexn3;
  radix_sort(indexn_input, indexnsize);

  int Last = 0;
  indexn3[indexn_input[0]] = 0;
  for (int i = 1; i < indexnsize; i++) {
    if (indexn_input[i] != indexn_input[Last]) {
      indexn_input[++Last] = indexn_input[i];
      indexn3[indexn_input[Last]] = Last;
    }
  }

  for (int i = 0; i < indexnsize; i++) {
    indexn_transform[i] = indexn3[indexn_transform[i]];
  }

  indexn_transform.erase(indexn_transform.begin() + Last,
                         indexn_transform.end());
  indexn_input.erase(indexn_input.begin() + Last, indexn_input.end());
  indexnsize = indexn_transform.size();
  /*
  *******************
  Windowing
  *******************
  */
  std::vector<complex_t> XW(np_channels * indexnsize);

  for (int i = 0; i < indexnsize; i++) {
    for (int j = 0; j < np_channels; j++) {
      XW[j + i * np_channels] =
          input_signal[indexn_input[i] + j] * chebwin_64[j];
    }
  }

  /*
  *******************
  First FFT
  *******************
  */
  std::vector<complex_t> XF1(np_channels * indexnsize);

  std::vector<complex_t> input_XFx1(np_channels);
  std::vector<complex_t> output_XFf1(np_channels);
  real_t(*ptr_input)[2] = reinterpret_cast<real_t(*)[2]>(input_XFx1.data());
  real_t(*ptr_output)[2] = reinterpret_cast<real_t(*)[2]>(output_XFf1.data());

  fftw_plan p;

  p = fftw_plan_dft_1d(np_channels, ptr_input, ptr_output, FFTW_FORWARD,
                       FFTW_MEASURE);
  for (int j = 0; j < indexnsize; j++) {
    for (int i = 0; i < np_channels; i++) {
      input_XFx1[i] = XW[i + j * np_channels];
    }

    fftw_execute(p);
    for (int i = 0; i < np_channels / 2; i++) {
      XF1[i + j * np_channels] =
          output_XFf1[i + np_channels / 2];  // with frequency shift
    }
    for (int i = np_channels / 2; i < np_channels; i++) {
      XF1[i + j * np_channels] =
          output_XFf1[i - np_channels / 2];  // with frequency shift
    }
  }
  fftw_destroy_plan(p);

  /*
  *******************************************
  Down Conversion + Conjugate Multiplication
  *******************************************
  */

  std::vector<complex_t> XM(indexnsize);
  std::vector<complex_t> Et;

  for (int k = 0; k < np_channels; k++) {
    real_t coeff_idx = real_t(k) / np_channels;
    Et.push_back(
        complex_t(cos(-2.0 * M_PI * coeff_idx), sin(-2.0 * M_PI * coeff_idx)));
  }

  /*
  ************************
  Second FFT
  ************************
  */
  /*prepare for sfft*/
  std::map<unsigned int, real_t> ans;
  SfftFuncParams sfft_func_parameters;
  sfft_func_parameters.size_of_n = size_of_n;
  sfft_func_parameters.np_channels = np_channels;
  sfft_func_parameters.kappa = kappa;
  sfft_func_parameters.b_est = B_est;
  sfft_func_parameters.b_thresh = B_thresh;
  sfft_func_parameters.b_loc = B_loc;
  sfft_func_parameters.parameters = parameters;
  sfft_func_parameters.index_trans = &indexn_transform;
  sfft_func_parameters.sigma_randn = &sigma_randn;
  sfft_func_parameters.myversion = myversion;
  sfft_func_parameters.filter_location = &filter_location;
  sfft_func_parameters.filter_estimation = &filter_estimation;
  sfft_func_parameters.size_of_comb_window = w_comb;

  for (int j = 0; j < np_channels; j++) {
    for (int m = 0; m < indexnsize; m++) {
      int myIndex = indexn_input[m];
      complex_t xc = conj(input_signal[myIndex + np_channels / 2]);
      XM[m] =
          Et[(((j + np_channels / 2) % np_channels) * myIndex) % np_channels] *
          XF1[j + m * np_channels] * xc;
    }

    std::map<unsigned int, real_t> an = SFFTfunc33(XM, sfft_func_parameters, j);
    ans.insert(an.begin(), an.end());
  }

  return ans;
}

std::map<unsigned int, real_t> single_S3CA(
    const std::vector<complex_t> &input_signal, int np_channels, int size_of_n,
    int seeds) {
  // Sparse FFT parameter library based on kappa = 50;
  // kappa is the expected sparse coefficient in the range of size_of_n;
  // kappa should be large or equal the expected number of coefficients in the
  // range;
  // If unkown the number of coefficient values, recommand to start from 50;
  constexpr int kDefaultKappa = 50;

  SfftVersion myversion;
  myversion.with_comb_filter = 1;
  myversion.seed = seeds;
  SfftParameters parameters =
      GetParameters(size_of_n, myversion.with_comb_filter);

  Filter filter_location =
      make_filter(size_of_n, parameters.loc_bcst, kDefaultKappa,
                  parameters.tolerance_loc, 1.1 * 1.2);
  Filter filter_estimation =
      make_filter(size_of_n, parameters.est_bcst, kDefaultKappa,
                  parameters.tolerance_est, 1.1 * 1.4);

  std::map<unsigned int, real_t> ans = sparsessca7(
      myversion, input_signal, np_channels, size_of_n, kDefaultKappa,
      parameters, filter_location, filter_estimation);

  return ans;
}

std::unique_ptr<std::vector<real_t>> single_SSCA(
    const std::vector<complex_t> &input_signal, int np_channels,
    int size_of_n) {
  std::unique_ptr<std::vector<real_t>> Sx =
      autossca_fftw(input_signal, np_channels, size_of_n);
  return Sx;
}

}  // namespace cyclostationary