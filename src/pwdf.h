#ifndef PWDF_H
#define PWDF_H

#include <complex>
#include <iostream>
#include <vector>

using std::complex;
using std::vector;

class PolynomialWDF4 {
   public:
    PolynomialWDF4(double d1 = -1.1151, double d2 = 0.4484, double d1m = 0.6667, double d2m = 1.0);
    vector<vector<complex<double>>> operator()(const vector<complex<double>>& x, double dt, double df);

   private:
    const double d1, d2, d1m, d2m;
    const double d[4] = {d1, d2, d1m, d2m};

    int p_min(int n, int n1, int n2);
    int p_max(int n, int n1, int n2);
    complex<double> interp(const vector<complex<double>>& x, double target);
};

#endif
