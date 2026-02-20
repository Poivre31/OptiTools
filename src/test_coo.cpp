#include <iostream>
#include "solver/solver.h"
#include "solver/coordinates.h"
#include <fmt/ranges.h>
#include <matplot/matplot.h>
#include <Eigen/Core>
#include <sys/resource.h>
#include <chrono>

using namespace std::chrono;

using namespace matplot;

double m = 1.;
double M = 7e22;
double G = 6.674e-11;

vec2 dvdt(double t, vec2 x)
{
    double r = x.norm();
    if (r == 0)
    {
        fmt::println("ERROR: division by 0 incoming, simulation invalidated");
        r = 1;
    }
    return -G * M * x / (r * r * r);
};

double potential(double t, vec2 x)
{
    double r = x.norm();
    if (r == 0)
    {
        fmt::println("ERROR: division by 0 incoming, simulation invalidated");
        r = 1;
    }
    return -G * M * m / r;
};

int main()
{
    auto t1 = high_resolution_clock::now();

    solver_degree_II<vec2> solver;
    double r0 = 2e6;
    double v = sqrt(G * M / r0);
    double T = 2 * pi * sqrt(r0 * r0 * r0 / (G * M));
    solver.set_initial_state(0, {r0, 0.}, {0., v});
    double tf = 20000000 * T;
    solver.set_timestep(T / 50);

    // plot(parse(solver.get_positions(), 0), parse(solver.get_positions(), 1));
    // show();
    // solver.solve_verlet(tf, dvdt);
    // auto E1 = solver.get_energy(m, potential);
    solver.solve_verlet(tf, dvdt);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);

    auto E2 = solver.get_energy(m, potential);
    // solver.solve_RK4(tf, dvdt);
    // auto E3 = solver.get_energy(m, potential);
    fmt::println("Exact energy : {}, computed: {}", -G * M * m / (2 * r0), E2[1]);
    // hold(on);
    // double E0 = E2[0];
    // E2 = apply_element_wise<double, double>(E2, [E0](double E)
    //                                         { return E - E0; });
    // // plot(solver.get_timeline(tf), E1)->display_name("Verlet");
    // // plot(stride(solver.get_timeline(tf), 100000), stride(E2, 100000))->display_name("Yoshida 4th");
    // // plot(solver.get_timeline(tf), E3)->display_name("RK4");
    // matplot::legend();
    // show();
}