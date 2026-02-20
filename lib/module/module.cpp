#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <solver/euler.h>

namespace nb = nanobind;

NB_MODULE(opti_tools, m)
{
    m.doc() = "Mon module de controle optimal";

    // m.def("solve_euler", &solve_euler, nb::arg("t0"), nb::arg("tf"), nb::arg("N"), nb::arg("x0"), nb::arg("v"), "Calcule x(tf) pour l'état initial x0 sous la contrainte dx/dt=v(t,x)");

    // nb::class_<solver>(m, "solver")
    //     .def(nb::init<size_t>())
    //     .def("solve", &solver::solve);
    m.def("add", [](int a, int b)
          { return a + b; });
    m.def("sub", [](int a, int b)
          { return a - b; });

    m.attr("__version__") = "dev";
}