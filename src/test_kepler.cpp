#include "solver/solver.h"
#include "solver/coordinates.h"
#include "celest/kepler_orbit.h"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <matplot/matplot.h>

constexpr double m = 1;
constexpr double M = 6e24;
constexpr double G = 6.67430e-11;

namespace plt = matplot;

xy a(double t, xy pos)
{
    double r = pos.norm();
    return -G * M * pos / (r * r * r);
}

int main()
{
    auto fig = plt::figure();
    fig->size(1000, 1000);
    plt::hold(plt::on);

    solver_degree_II<xy> solver;
    xy r0 = {6371000., 2371000.};
    xy v0 = {.8 * sqrt(G * M / r0.norm()), .9 * sqrt(G * M / r0.norm())};
    // auto marker1 = plt::plot({r0[0]}, {r0[1]}, ".");
    // marker1->marker_color({0.8, 0.4, 0.4});
    // marker1->display_name("r(t0)");

    orbit orbit(M, m, G);
    orbit.set_from_current_state(r0, v0);

    solver.set_initial_state(0, r0, v0);
    solver.set_timestep(orbit.get_period() / 1000);
    std::vector<xy> XY;

    double tf = 1 * orbit.get_period();
    double pad = .1 * orbit.get_semi_major_axis();
    plt::xlim({-pad - orbit.get_r_max(), pad + orbit.get_r_min()});
    plt::ylim({-pad - orbit.get_semi_major_axis(), pad + orbit.get_semi_major_axis()});
    auto marker2 = plt::plot({0.}, {0.}, ".");
    marker2->marker_color({0.1, 0.6, 0.5});
    marker2->display_name("Earth");

    solver.solve_RK4(tf, a);
    XY = solver.get_positions();
    plt::plot(parse(XY, 0), parse(XY, 1), "--")->display_name("RK4");

    auto exact_trajectory = orbit.get_trajectory(1000);
    auto plot = plt::plot(parse(exact_trajectory, 0), parse(exact_trajectory, 1));
    plot->display_name("Exact trajectory");
    plot->color("blue");

    plt::legend();
    plt::hold(plt::off);
    plt::show();
}