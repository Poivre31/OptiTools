import opti_tools


def v(t, x):
    return 0.1 * x


solver = opti_tools.solver(4)
print(solver.solve(3))
