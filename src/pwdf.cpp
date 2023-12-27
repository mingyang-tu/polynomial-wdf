#include "pwdf.h"

inline void fft1d(vcd1d& x, vcd1d& y, int N) {
    fftw_complex* in = reinterpret_cast<fftw_complex*>(x.data());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(y.data());
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

PolynomialWDF4::PolynomialWDF4(double d1, double d2, double d1m, double d2m)
    : d1(d1), d2(d2), d1m(d1m), d2m(d2m) {
    double eps = 1e-3;
    double rule1 = d1 + d2 + d1m + d2m - 1;
    double rule2 = pow(d1, 2) + pow(d2, 2) - pow(d1m, 2) - pow(d2m, 2);
    double rule3 = pow(d1, 3) + pow(d2, 3) + pow(d1m, 3) + pow(d2m, 3);

    if (abs(rule1) > eps)
        cout << "Warning: d_1 + d_2 + d_-1 + d_-2 = " << rule1 << " (should be 1)" << endl;
    if (abs(rule2) > eps)
        cout << "Warning: d_1^2 + d_2^2 - d_-1^2 - d_-2^2 = " << rule2 << " (should be 0)" << endl;
    if (abs(rule3) > eps)
        cout << "Warning: d_1^3 + d_2^3 + d_-1^3 + d_-2^3 = " << rule3 << " (should be 0)" << endl;
}

// x: input signal, t: time vector, f: frequency vector
// dt: time step, df: frequency step
// return: 2D vector of complex numbers
vcd2d PolynomialWDF4::operator()(const vcd1d& x, const vector<double>& t, const vector<double>& f,
                                 double dt, double df) {
    if (t.size() != x.size())
        throw invalid_argument("t and x must have the same size");

    const int T = t.size();
    const int F = f.size();
    const int N = int(1 / (dt * df));
    const int n1 = 0, n2 = x.size() - 1;

    vcd1d x_conj(T);
    for (int i = 0; i < T; i++)
        x_conj[i] = conj(x[i]);

    int f_start = (int(f[0] / df) % N + N) % N;
    vector<int> m(F);
    for (int i = 0; i < F; i++)
        m[i] = (i + f_start) % N;

    vcd2d output(F, vcd1d(T));
    complex<double> bias = complex<double>(0, 1) * (2.0 * M_PI / N);

    for (int n = 0; n < T; n++) {
        vcd1d kernel(N, 0);
        int p_min = this->p_min(n, n1, n2), p_max = this->p_max(n, n1, n2);
        for (int p = p_min; p <= p_max; p++) {
            kernel[p - p_min] = interp(x, n + d1 * p) * interp(x_conj, n - d1m * p) *
                                interp(x, n + d2 * p) * interp(x_conj, n - d2m * p);
        }
        vcd1d dft(N);
        fft1d(kernel, dft, N);
        complex<double> bias_p = bias * (double)p_min;
        for (int i = 0; i < F; i++) {
            output[i][n] = dt * exp(bias_p * double(m[i])) * dft[m[i]];
        }
    }

    return output;
}

int PolynomialWDF4::p_min(int n, int n1, int n2) {
    double p_min = n1 - n2;
    for (int i = 0; i < 4; i++) {
        if (d[i] > 0)
            p_min = max((n1 - n) / d[i], p_min);
        else if (d[i] < 0)
            p_min = max((n2 - n) / d[i], p_min);
    }
    return ceil(p_min);
}

int PolynomialWDF4::p_max(int n, int n1, int n2) {
    double p_max = n2 - n1;
    for (int i = 0; i < 4; i++) {
        if (d[i] > 0)
            p_max = min((n2 - n) / d[i], p_max);
        else if (d[i] < 0)
            p_max = min((n1 - n) / d[i], p_max);
    }
    return floor(p_max);
}

complex<double> PolynomialWDF4::interp(const vcd1d& x, double target) {
    int x_floor = floor(target);
    int x_ceil = ceil(target);
    return x[x_floor] + (x[x_ceil] - x[x_floor]) * (target - x_floor);
}
