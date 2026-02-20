#pragma once
#include <vector>
#include <functional>
#include "solver/coordinates.h"

template <typename T>
std::vector<T> midpoint(
    double t0, double tf, size_t N, T x0, std::function<T(double t, T x)> dxdt)
{
    double dt = (tf - t0) / N; ///< We consider N+1 steps from t0 to tf included so dt = T/N
    T x = x0;

    std::vector<T> positions(N + 1);
    positions[0] = x0;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        auto vi = dxdt(t, x);
        x += dt * dxdt(t + dt / 2, x + dt / 2 * vi);

        positions[i] = x;

        t += dt;
    }
    return positions;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> midpoint(
    double t0, double tf, size_t N, T x0, T v0, std::function<T(double t, T x, T v)> a)
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
        T x_i = x + dt / 2 * v;
        T v_i = v + dt / 2 * a(t, x, v);
        x += dt * v_i;
        v += dt * a(t + dt / 2, x_i, v_i);

        positions[i] = x;
        velocities[i] = v;

        t += dt;
    }
    return std::make_tuple(positions, velocities);
}