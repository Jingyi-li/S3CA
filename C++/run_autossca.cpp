#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <map>
#include <cassert>
#include <algorithm>
#include <cstring>

#include "SSCA.h"
#include "types.h"
#include "fftw.h"
#include "computefourier.h"
#include "utils.h"
#include "autossca.h"

using namespace std;
using namespace cyclostationary;
// Use the autossca library compate the L1 norm with different K for FFTW and Sparse FFT
// FFTW will pick first lagest K value like sparse FFT
/*
  x: the signal
  n: the length of x

  lobefrac_loc:   during location, the main lobe of the filter has half
                  width n*lobefrac_loc
  tolerance_loc:  the linf norm of the residuals of the filter
  b_loc:          the number of adjacent filters to add

  B_loc:          number of samples in subsampling during location
  B_thresh:       number of samples considered "heavy"
  loops_loc:      number of location loops
  loops_thresh:   number of times a coordinate must seem heavy.

  *_est:          ditto as above, but for estimation.

  repetitions:    repeat the experiment this many times for timing // no repetitions
  LARGE_FREQ:     locations of largest coefficients.
  k:              number of HHs, used for evaluation only.
  x_f:            true fft.
 */

int main()
{
   int data_length = 400000, size_of_n = 65536, np_channels = 64;
  std::vector<complex_t> input_signal;
  FILE* fp;
  fp = fopen("verify_data.dat", "r");
  if (fp == 0) {printf("Error: no matched file.\n");assert(0);}

  real_t value_i, value_q;
  int mylength = size_of_n > data_length ? data_length : size_of_n;
  for (int ii = 0; ii < mylength; ii++) {
    if (std::fscanf(fp, "%f\t%f\n", &value_i, &value_q) != 0) {
      input_signal.push_back(complex_t(value_i, value_q));
    } else {
      printf("Error: file length not match.\n");assert(0);
    }
  }
  for (int ii = mylength; ii < size_of_n + np_channels; ii++) {
    input_signal.push_back(complex_t(0, 0));
  }
  fclose(fp);
  assert(input_signal.size()== (size_of_n + np_channels));

  // Run SSCA
  std::unique_ptr<std::vector<real_t>> Sx =
      single_SSCA(input_signal, np_channels, size_of_n);
  std::vector<real_t>& Sx_vec = *Sx.get();
  assert(Sx_vec.size() == (np_channels * size_of_n * 2));
  assert(fabs(Sx_vec[3774882] - 2.3541) < 0.0001);
  assert(fabs(Sx_vec[4194344] - 1.4982) < 0.0001);

  //Run S3CA
  std::map<unsigned int, real_t> Sx_S3CA =
      single_S3CA(input_signal, np_channels, size_of_n, 23);

  assert(fabs(Sx_S3CA[3355487] - 1.170636) < 0.0001);
  assert(fabs(Sx_S3CA[3774882] - 2.440637) < 0.0001);
  assert(fabs(Sx_S3CA[4194343] - 3.192273 ) < 0.0001);
  assert(fabs(Sx_S3CA[4194338] - 8.362414) < 0.0001);
  printf("DONE\n");
    return 0;
}