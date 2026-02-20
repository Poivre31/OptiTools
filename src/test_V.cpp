#include <iostream>
#include "solver/solver.h"
#include <fmt/ranges.h>
#include <matplot/matplot.h>
#include "math/vec2d.h"
#include "solver/coordinates.h"

using namespace matplot;

int main()
{
    double q = 1;
    double ε0 = 1e-3;
    double π = M_PI;
    double a = 0.01;
    vec2d P{a, 0.};
    vec2d N{-a, 0.};

    auto E = [P, N, q, ε0, π](double, vec2d pos)
    {
        return q / (4 * π * ε0) * ((pos - P) / pow((pos - P).norm(), 3) - (pos - N) / abs(pow((pos - N).norm(), 3)));
    };

    solver_degree_1<vec2d> solver;
    vec2d x0 = {1., 1.};
    solver.set_initial_state(0, x0);
    solver.set_timestep(0.01);
    solver.solve_RK4(100, E);
    auto pos = solver.get_positions();
    std::vector<double> X(pos.size());
    std::vector<double> Y(pos.size());
    for (size_t i = 0; i < pos.size(); i++)
    {
        X[i] = pos[i].x;
        Y[i] = pos[i].y;
    }

    double C = pow(x0.norm(), 3) / (x0.y * x0.y);
    // hold(on);
    // plot(X, Y);
    // fimplicit([C](double x, double y)
    //           { return pow(x * x + y * y, 1.5) - C * y * y; });
    // xlim({-5 * x0.x, 5 * x0.x});
    // ylim({-5 * x0.y, 5 * x0.y});
    // show();

    vec3d aaa;
    aaa.x = 3;
    aaa.data[1] = 7;
    aaa.z = 2;
    std::cout << aaa.data << std::endl;
}