#ifndef __CYCLOSTATIONARY_AUTOSSCA_H__
#define __CYCLOSTATIONARY_AUTOSSCA_H__

#include <random>

#include "computefourier.h"
#include "filters.h"
#include "parameters.h"
#include "types.h"

namespace cyclostationary {

/*Compute the dolphchebyshev filter with time domain and frequency domain. Time
 * domain save the support size, and frequency domain save the full size*/
Filter make_filter(int size_of_n, double Bcst, int kappa, double tolerance,
                   double factor);

/*This function compute the sparse fft based on given parameters*/
std::map<unsigned int, real_t> SFFTfunc33(
    std::vector<complex_t> &input_signal,
    const SfftFuncParams sfft_func_parameters, int indx_channel);

/*Compute the SSCA based on fftw library*/
std::unique_ptr<std::vector<real_t>> autossca_fftw(
    const std::vector<complex_t> &input_signal, int np_channels, int size_of_n);
/*Generate random offset for comb loops and random co prime sigma for location
 * and estimation loops*/
std::vector<int> GenerateSigmaRandom(std::mt19937 &generator, int size_of_n,
                                     int size_of_comb, int comb_loops,
                                     int sumloops);
/*****************sparsessca7********************
//1. Base on the sparse FFT compute the sparse
//   window in first section
//2. Consider the filter function as an input
//3. Only output non-zero values with location info
//4. Pass the location information of subset X_g to sparse fft to have the
//same location information of N X_g
***************************************************/
std::map<unsigned int, real_t> sparsessca7(
    const SfftVersion &myversion, const std::vector<complex_t> &input_signal,
    int np_channels, int size_of_n, int kappa, const SfftParameters &parameters,
    const Filter &filter_location, const Filter &filter_estimation);
/*package sparsessca in single version*/
std::map<unsigned int, real_t> single_S3CA(
    const std::vector<complex_t> &input_signal, int np_channels, int size_of_n,
    int seeds);
/*package ssca in single version*/
std::unique_ptr<std::vector<real_t>> single_SSCA(
    const std::vector<complex_t> &input_signal, int np_channels, int size_of_n);

}  // namespace cyclostationary
#endif  // __CYCLOSTATIONARY__AUTOSSCA_H__
