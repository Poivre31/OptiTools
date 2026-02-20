#include "euler.h"
#include <math/vec3d.h>

double solve_euler(double t0, double tf, size_t N, double x0, std::function<double(double, double)> v)
{
    double dt = (tf - t0) / N;
    double x = x0;

    double t = dt;
    for (size_t i = 1; i <= N; i++)
    {
        x += dt * v(t, x);
        t += dt;
    }
    return x;
}