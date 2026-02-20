#pragma once
#include <vector>
#include <functional>
#include "solver/coordinates.h"
#include <iostream>

template <typename T>
std::vector<T> euler_explicit(
    double t0, double tf, size_t N, T x0, std::function<T(double t, T x)> dxdt)
{
    double dt = (tf - t0) / N; ///< We consider N+1 steps from t0 to tf included so dt = T/N
    T x = x0;

    std::vector<T> positions(N + 1);
    positions[0] = x0;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        x += dt * dxdt(t, x);

        positions[i] = x;

        t += dt;
    }
    return positions;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> euler_explicit(
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
        T a_n = a(t, x, v);
        x += dt * v;
        v += dt * a_n;

        positions[i] = x;
        velocities[i] = v;

        t += dt;
    }
    return std::make_tuple(positions, velocities);
}