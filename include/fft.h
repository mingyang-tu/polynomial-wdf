#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>

#include "fftw3.h"

using std::complex;
using std::vector;

void fft1d(vector<complex<double>>& x, vector<complex<double>>& y) {
    const int N = x.size();
    fftw_complex* in = reinterpret_cast<fftw_complex*>(x.data());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(y.data());
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

#endif
