import numpy as np
from math import ceil, floor
import matplotlib.pyplot as plt
import time


class PolynomialWDF4:
    def __init__(self, d1=-1.1151, d2=0.4484, d1m=0.6667, d2m=1.0):
        eps = 1e-3
        rule1 = d1 + d2 + d1m + d2m - 1
        rule2 = d1**2 + d2**2 - d1m**2 - d2m**2
        rule3 = d1**3 + d2**3 + d1m**3 + d2m**3
        if abs(rule1) > eps:
            print(f"Warning: d_1 + d_2 + d_-1 + d_-2 = {rule1} (should be 1)")
        if abs(rule2) > eps:
            print(f"Warning: d_1^2 + d_2^2 - d_-1^2 - d_-2^2 = {rule2} (should be 0)")
        if abs(rule3) > eps:
            print(f"Warning: d_1^3 + d_2^3 + d_-1^3 + d_-2^3 = {rule3} (should be 0)")
        self.d1 = d1
        self.d2 = d2
        self.d1m = d1m
        self.d2m = d2m
        self.d_all = [d1, d2, -d1m, -d2m]

    def __call__(self, x, t, f, dt, df):
        assert len(x) == len(t)
        T = len(t)
        F = len(f)
        N = int(1 / (dt * df))
        n1, n2 = 0, len(x) - 1
        f1 = int(f[0] / df)
        f2 = f1 + F - 1

        x_conj = np.conj(x)
        f_index = np.arange(f1, f2 + 1) % N
        m = np.arange(0, N, dtype=np.float64)[f_index]
        output = np.zeros((F, T), dtype=np.complex128)

        for n in range(T):
            kernel = np.zeros(N, dtype=np.complex128)
            p_min, p_max = self.p_range(n, n1, n2)
            p_indexes = np.arange(p_min, p_max + 1)
            kernel[: p_indexes.shape[0]] = (
                self.interp(x, n + self.d1 * p_indexes)
                * self.interp(x_conj, n - self.d1m * p_indexes)
                * self.interp(x, n + self.d2 * p_indexes)
                * self.interp(x_conj, n - self.d2m * p_indexes)
            )
            dft = np.fft.fft(kernel)
            output[:, n] = dt * np.exp(1j * 2 * np.pi * m * p_min / N) * dft[f_index]

        return output

    def p_range(self, n, n1, n2):
        p_max = n2 - n1
        p_min = -p_max
        for d in self.d_all:
            if d > 0:
                p_min = max((n1 - n) / d, p_min)
                p_max = min((n2 - n) / d, p_max)
            elif d < 0:
                p_min = max((n2 - n) / d, p_min)
                p_max = min((n1 - n) / d, p_max)
        return ceil(p_min), floor(p_max)

    def interp(self, x, target):
        x_floor = np.floor(target).astype(int)
        x_ceil = np.ceil(target).astype(int)
        return x[x_floor] + (x[x_ceil] - x[x_floor]) * (target - x_floor)


def visualize(wigner, t, f):
    t_min, t_max = t[0], t[-1]
    f_min, f_max = f[0], f[-1]

    plt.figure()
    plt.imshow(
        np.abs(wigner),
        cmap="gray",
        origin="lower",
        aspect="auto",
        extent=(t_min, t_max, f_min, f_max),
    )
    plt.title("Polynomial Wigner Distribution Function")
    plt.xlabel("Time")
    plt.ylabel("Frequency")
    plt.show()


def main():
    dt = 0.0125
    df = 0.0125

    t = np.arange(0, 10, dt)
    f = np.arange(-5, 15, df)
    x = np.exp(1j * (t - 5) ** 3)

    start = time.time()
    pwdf = PolynomialWDF4()(x, t, f, dt, df)
    end = time.time()
    print(f"Elapsed time: {end - start} s")

    visualize(pwdf, t, f)


if __name__ == "__main__":
    main()
