#pragma once
#include <vector>
#include <functional>
#include "solver/coordinates.h"
#include <iostream>

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> yoshida_4th(
    double t0, double tf, size_t N, T x0, T v0, std::function<T(double t, T x)> a)
{
    double dt = (tf - t0) / N; ///< We consider N+1 steps from t0 to tf included so dt = T/N
    T x = x0;
    T v = v0;

    std::vector<T> positions(N + 1);
    positions[0] = x0;
    std::vector<T> velocities(N + 1);
    velocities[0] = v0;

    double ω0 = -cbrt(2) / (2 - cbrt(2));
    double ω1 = 1 / (2 - cbrt(2));
    double c1 = ω1 / 2;
    double c2 = (ω0 + ω1) / 2;
    double c3 = (ω0 + ω1) / 2;
    double c4 = ω1 / 2;
    double d1 = ω1;
    double d2 = ω0;
    double d3 = ω1;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        // if (i % 10000000 == 0)
        //     std::cout << 100. * i / N << "\n";
        T x1 = x + c1 * v * dt;
        T v1 = v + d1 * a(t + c1 * dt, x1) * dt;
        T x2 = x1 + c2 * v1 * dt;
        T v2 = v1 + d2 * a(t + (c1 + c2) * dt, x2) * dt;
        T x3 = x2 + c3 * v2 * dt;
        T v3 = v2 + d3 * a(t + (c1 + c2 + c3) * dt, x3) * dt;

        x = x3 + c4 * v3 * dt;
        v = v3;

        positions[i] = x;
        velocities[i] = v;

        t += dt;
    }

    return std::make_tuple(positions, velocities);
}