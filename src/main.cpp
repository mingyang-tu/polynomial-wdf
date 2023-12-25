#include "matplotlibcpp.h"
#include "pwdf.h"
#include "utils.h"

#define JJ complex<double>(0, 1)

namespace plt = matplotlibcpp;

using namespace std;

// visualize the result
void visualize(vector<vector<complex<double>>> result, double t_min, double t_max, double f_min,
               double f_max, double df) {
    int f_idx_start = int(f_min / df), f_idx_end = int(f_max / df);
    int result_size = result.size();
    int nrows = f_idx_end - f_idx_start, ncols = result[0].size();

    vector<float> image(ncols * nrows);
    for (int i = 0; i < nrows; i++) {
        int i_result = ((i + f_idx_start) % result_size + result_size) % result_size;
        for (int j = 0; j < ncols; j++) {
            image[ncols * i + j] = abs(result[i_result][j]);
        }
    }

    const float* image_ptr = &(image[0]);

    plt::title("Polynomial Wigner Distribution Function");
    plt::imshow(image_ptr, nrows, ncols, 1, {{"cmap", "gray"}, {"origin", "lower"}, {"aspect", "auto"}});
    plt::xlabel("Time");
    plt::ylabel("Frequency");

    int num_ticks_x = 5, num_ticks_y = 5;
    vector<int> xticks = get_ticks(0, ncols, num_ticks_x), yticks = get_ticks(0, nrows, num_ticks_y);
    vector<string> xticklabels = get_ticklabels(t_min, t_max, num_ticks_x),
                   yticklabels = get_ticklabels(f_min, f_max, num_ticks_y);
    plt::xticks(xticks, xticklabels);
    plt::yticks(yticks, yticklabels);
    plt::show();
}

int main() {
    double dt = 0.0125, df = 0.0125;
    double t_start = 0, t_end = 10;
    const int T = static_cast<int>((t_end - t_start) / dt);

    vector<complex<double>> x(T);
    for (int i = 0; i < T; i++) {
        double t = t_start + i * dt;
        x[i] = exp(JJ * pow(t - 5, 3));
    }

    PolynomialWDF4 pwdf;

    Timer clock;
    clock.tick();

    vector<vector<complex<double>>> result = pwdf(x, dt, df);

    clock.tock();
    cout << "Elapsed time: " << clock.duration() << " s" << endl;

    visualize(result, t_start, t_end, -5, 15, df);

    return 0;
}