#pragma once
#include <vector>
#include <functional>
#include "solver/coordinates.h"

template <typename T>
std::vector<T> RK4_explicit(
    double t0, double tf, size_t N, T x0, std::function<T(double t, T x)> dxdt)
{
    double dt = (tf - t0) / N; ///< We consider N+1 steps from t0 to tf included so dt = T/N
    T x = x0;

    std::vector<T> positions(N + 1);
    positions[0] = x0;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        auto k1 = dxdt(t, x);
        auto k2 = dxdt(t + dt / 2, x + k1 * dt / 2);
        auto k3 = dxdt(t + dt / 2, x + k2 * dt / 2);
        auto k4 = dxdt(t + dt, x + k3 * dt);
        x += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

        positions[i] = x;

        t += dt;
    }
    return positions;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> RK4_explicit(
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
        auto v1 = v;
        auto x1 = x;
        auto k1 = a(t, x, v);

        auto v2 = v1 + dt / 2 * k1;
        auto x2 = x1 + dt / 2 * v1;
        auto k2 = a(t + dt / 2, x2, v2);

        auto v3 = v1 + dt / 2 * k2;
        auto x3 = x1 + dt / 2 * v2;
        auto k3 = a(t + dt / 2, x3, v3);

        auto v4 = v1 + dt * k3;
        auto x4 = x1 + dt * v3;
        auto k4 = a(t + dt, x4, v4);

        v += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        x += dt / 6 * (v1 + 2 * v2 + 2 * v3 + v4);

        positions[i] = x;
        velocities[i] = v;

        t += dt;
    }
    return std::make_tuple(positions, velocities);
}