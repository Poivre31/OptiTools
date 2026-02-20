#include <iostream>
#include "solver/solver.h"
#include <fmt/ranges.h>
#include <matplot/matplot.h>

using namespace matplot;

int main()
{
  solver_degree_1<double> solver;
  solver.set_initial_state(0, 3);
  solver.set_timestep(0.1);
  auto dxdt = [](double t, double x)
  { return 3 * x; };

  solver.solve_euler(1, dxdt);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver.get_positions(), ", "));

  solver.solve_midpoint(1, dxdt);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver.get_positions(), ", "));

  solver.solve_RK4(1, dxdt);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver.get_positions(), ", "));

  fmt::println("Analytical solution: {:.3g}", 3 * exp(3));

  solver_degree_2<double> solver_2;
  solver_2.set_initial_state(0, 3, 0);
  solver_2.set_timestep(0.1);
  double ω0 = 2;
  double τ = 3;
  double tf = 10;
  auto a = [ω0, τ](double t, double x, double v)
  {
    return -ω0 * ω0 * x - v / τ;
  };

  solver_2.solve_euler(tf, a);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver_2.get_positions(), ","));

  solver_2.solve_midpoint(tf, a);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver_2.get_positions(), ","));

  solver_2.solve_RK4(tf, a);
  fmt::println("Trajectory: {:.3g}", fmt::join(solver_2.get_positions(), ","));

  double ω = sqrt(ω0 * ω0 - 1 / (4 * τ * τ));
  fmt::println("Analytical solution: {:.3g}", exp(-tf / (2 * τ)) * (3 * cos(ω * tf) + 3 / (2 * τ * ω) * sin(ω * tf)));

  hold(on);
  plot(solver_2.get_timeline(tf), solver_2.get_positions(), solver_2.get_timeline(tf), solver_2.get_velocities());
  show();
}