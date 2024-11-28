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

#include "parameters.h"

namespace cyclostationary {
SfftParameters GetParameters(const int size_of_fft,
                             const bool with_comb_filter) {
  SfftParameters parameters;
  if (with_comb_filter) {
    switch (size_of_fft) {
      case 4096:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 2;
        parameters.comb_cst = 24;
        parameters.comb_loops = 4;
        parameters.est_loops = 8;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 8192:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 2;
        parameters.comb_cst = 32;
        parameters.comb_loops = 8;
        parameters.est_loops = 16;
        parameters.loc_loops = 7;
        parameters.threshold_loops = 6;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 16384:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 4;
        parameters.comb_cst = 32;
        parameters.comb_loops = 8;
        parameters.est_loops = 10;
        parameters.loc_loops = 6;
        parameters.threshold_loops = 5;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 32768:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 2;
        parameters.comb_cst = 64;
        parameters.comb_loops = 4;
        parameters.est_loops = 8;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 65536:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 2;
        parameters.comb_cst = 128;
        parameters.comb_loops = 6;
        parameters.est_loops = 10;
        parameters.loc_loops = 4;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 131072:
        parameters.loc_bcst = 1;
        parameters.est_bcst = 1;
        parameters.comb_cst = 8;
        parameters.comb_loops = 2;
        parameters.est_loops = 12;
        parameters.loc_loops = 4;
        parameters.threshold_loops = 3;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 262144:
        parameters.loc_bcst = 1;
        parameters.est_bcst = 1;
        parameters.comb_cst = 8;
        parameters.comb_loops = 2;
        parameters.est_loops = 14;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 524288:
        // loc_bcst =0.5; est_bcst =0.5; comb_cst =  8; comb_loops =1; est_loops
        // =10; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6;
        // tolerance_est =1e-8;
        parameters.loc_bcst = 1;
        parameters.est_bcst = 1;
        parameters.comb_cst = 8;
        parameters.comb_loops = 2;
        parameters.est_loops = 13;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 3;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 1048576:
        // loc_bcst =0.5; est_bcst =0.5; comb_cst =  8; comb_loops =2; est_loops
        // =12; loc_loops =4; threshold_loops =2; tolerance_loc =1e-6;
        // tolerance_est =1e-8;
        parameters.loc_bcst = 0.5;
        parameters.est_bcst = 0.5;
        parameters.comb_cst = 8;
        parameters.comb_loops = 1;
        parameters.est_loops = 10;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 1;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 2097152:
        parameters.loc_bcst = 0.5;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 8;
        parameters.comb_loops = 1;
        parameters.est_loops = 10;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 4194304:
        parameters.loc_bcst = 0.5;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 8;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 8388608:
        parameters.loc_bcst = 0.5;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 8;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 16777216:
        parameters.loc_bcst = 0.5;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 16;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
    }
  } else {
    switch (size_of_fft) {
      case 8192:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 16;
        parameters.loc_loops = 7;
        parameters.threshold_loops = 6;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 16384:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 4;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 10;
        parameters.loc_loops = 6;
        parameters.threshold_loops = 5;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 32768:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 65536:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-8;
        parameters.tolerance_est = 1e-8;
        break;
      case 131072:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 1;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 10;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 262144:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 0.5;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 14;
        parameters.loc_loops = 4;
        parameters.threshold_loops = 3;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 524288:
        parameters.loc_bcst = 1;
        parameters.est_bcst = 0.5;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 12;
        parameters.loc_loops = 5;
        parameters.threshold_loops = 4;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 1048576:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 0.5;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 12;
        parameters.loc_loops = 4;
        parameters.threshold_loops = 3;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 2097152:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 15;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 4194304:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 10;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 8388608:
        parameters.loc_bcst = 2;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
      case 16777216:
        parameters.loc_bcst = 4;
        parameters.est_bcst = 0.2;
        parameters.comb_cst = 1;
        parameters.comb_loops = 1;
        parameters.est_loops = 8;
        parameters.loc_loops = 3;
        parameters.threshold_loops = 2;
        parameters.tolerance_loc = 1e-6;
        parameters.tolerance_est = 1e-8;
        break;
    }
  }

  return parameters;
}

}  // namespace cyclostationary