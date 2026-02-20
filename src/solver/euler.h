#pragma once
#include <math/vec3d.h>
#include <functional>
#include <iostream>
#include <fmt/format.h>

double solve_euler(double t0, double tf, size_t N, double x0, std::function<double(double, double)> v);

class solver
{
public:
    solver(size_t degree) { set_degree(degree); };
    solver(size_t degree, double t0, double x0, double dt) : _x0(x0), _t0(t0)
    {
        set_degree(degree);
        set_stepsize(dt);
    };

    void set_degree(size_t n)
    {
        if (n < 1 || n > 2)
        {
            fmt::println("ERROR: Only degree 1 and 2 differential equations supported");
            return;
        }
        _degree = n;
    }

    void set_initial_conditions(double t0, double x0)
    {
        _t0 = t0;
        _x0 = x0;
    }

    void set_stepsize(double dt)
    {
        if (dt <= 0)
        {
        std:
            fmt::println("ERROR: Stepsize must be greater than 0");
            return;
        }
        _dt = dt;
    }

    double solve(double tf)
    {
        fmt::println("x0 = {0}, dt = {1}", _x0, _dt);
        return tf;
    }

private:
    size_t _degree = 1;
    double _t0 = 0;
    double _x0 = 0;
    double _dt = 1;
};