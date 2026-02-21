#include <solver/coordinates.h>
#include <solver/solver.h>
#include <matplot/matplot.h>
#include <functional>
#include <chrono>
#include <iostream>

using namespace std::chrono;

namespace plt = matplot;

using state = coordinates<5>;

namespace c
{
    enum
    {
        x,
        z,
        vx,
        vz,
        alpha,
    };
}

double M = 3e6;
double tf = 1000.;
double dt = 1;
double A = 1000 * 16.5;
double k = 2500 * A;
double g = 9.81;
double γ = 1.0e-5 * M;
double H = 8000.;

state dsdt(double t, state state_i)
{

    state state_change;
    double drag = γ * exp(-state_i[c::z] / H) / M;
    double vx = state_i[c::vx];
    double vz = state_i[c::vz];
    double v = sqrt(vx * vx + vz * vz);
    double α = state_i[c::alpha];

    state_change[c::x] = vx;
    state_change[c::z] = vz;
    state_change[c::vx] = k / M * sin(α) - drag * vx * v;
    state_change[c::vz] = k / M * cos(α) - drag * vz * v - g;
    if (v != 0)
    {
        state_change[c::alpha] = drag * (sin(2 * α) * (vx * vx - vz * vz) + cos(2 * α) * vx * vz) / (2. * v);
    }
    else
    {
        state_change[c::alpha] = 0;
    }

    return state_change;
}

std::vector<state> propagate(solver_degree_I<state> &solver, double α0)
{
    state X0 = {0., 0., 0., 0., α0};

    solver.set_initial_state(0, X0);
    solver.set_timestep(dt);
    solver.solve_RK4(tf, dsdt);
    auto result = solver.get_positions();
    size_t i_max = get_maximum_index(parse(result, c::z));
    return truncate_vector(result, i_max);
}

double gradient(solver_degree_I<state> &solver, double target, double α0, double dα)
{
    auto result = propagate(solver, α0 - 2 * dα);
    double h_m2 = result.back()[c::z];
    double L2_m2 = (h_m2 / target - 1) * (h_m2 / target - 1);

    result = propagate(solver, α0 - dα);
    double h_m1 = result.back()[c::z];
    double L2_m1 = (h_m1 / target - 1) * (h_m1 / target - 1);

    result = propagate(solver, α0 + dα);
    double h_p1 = result.back()[c::z];
    double L2_p1 = (h_p1 / target - 1) * (h_p1 / target - 1);

    result = propagate(solver, α0 + 2 * dα);
    double h_p2 = result.back()[c::z];
    double L2_p2 = (h_p2 / target - 1) * (h_p2 / target - 1);

    // fmt::println("{} {} {} {}", h_m2, h_m1, h_p1, h_p2);

    return (1. / 12 * L2_m2 - 2. / 3 * L2_m1 + 2. / 3 * L2_p1 - 1. / 12 * L2_p2) / dα;
}

int main()
{
    auto t1 = high_resolution_clock::now();
    double α0 = .9;

    solver_degree_I<state> solver;
    double dα = 0.000001;
    double target = 50000;

    double alpha = α0;
    double grad = gradient(solver, target, α0, dα);
    double n_alpha = alpha - grad * 1e-10;
    // fmt::println("alpha: {:.10f}, grad: {:.2g}", n_alpha, grad);
    for (size_t i = 0; i < 1000; i++)
    {
        double n_grad = gradient(solver, target, n_alpha, dα);
        if (n_grad == grad)
            break;
        double η = 0.1 * abs((n_alpha - alpha) / (n_grad - grad));
        grad = n_grad;
        alpha = n_alpha;
        n_alpha -= grad * η;
        double z = propagate(solver, n_alpha).back()[c::z];
        // fmt::println("I: {}, alpha: {:.10f}, eta: {:.2g}, Z: {:.6g}, grad: {:.2g}", i, n_alpha, η, z, grad);
        if (abs(z - target) < 0.5 || std::isnan(α0) || std::isnan(η))
            break;
    }

    auto result = propagate(solver, alpha);
    auto time = truncate_vector(solver.get_timeline(tf), result.size() + 1);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);
    fmt::println("Final height: {:.1f}m, velocity: {:.1f}m/s, time: {:.1f}s, conso: {:.1f}T", parse(result, c::z).back(), parse(result, c::vx).back(), time.back(), time.back() * A / 1000);

    if (false)
    {
        auto fig = plt::figure();
        fig->size(1920, 1080);
        plt::subplot(2, 2, 0);
        plt::plot(parse(result, c::x), parse(result, c::z));
        plt::title("Trajectory z(x)");

        plt::subplot(2, 2, 1);
        plt::plot(time, parse(result, c::alpha));
        plt::title("Angle α(t)");

        plt::subplot(2, 2, 2);
        plt::plot(time, parse(result, c::z));
        plt::title("Alitutde z(t)");

        plt::subplot(2, 2, 3);
        plt::plot(time, parse(result, c::vx));
        plt::title("Ortho velocity vx(t)");

        plt::show();
    }
}
