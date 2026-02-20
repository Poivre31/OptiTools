#pragma once
#include <vector>
#include <functional>
#include "solver/coordinates.h"
#include <iostream>

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> verlet_velocity(
    double t0, double tf, size_t N, T x0, T v0, std::function<T(double t, T x)> a)
{
    double dt = (tf - t0) / N; ///< We consider N+1 steps from t0 to tf included so dt = T/N
    T x = x0;
    T v = v0;

    std::vector<T> positions(N + 1);
    positions[0] = x0;
    std::vector<T> velocities(N + 1);
    velocities[0] = v0;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        T v_half = v + dt / 2 * a(t, x);
        x += dt * v_half;
        v = v_half + dt / 2 * a(t, x);

        positions[i] = x;
        velocities[i] = v;

        t += dt;
    }
    return std::make_tuple(positions, velocities);
}