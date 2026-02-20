import opti_tools

print(opti_tools.add(1))


def v(t, x):
    return 0.1 * x


def square(x):
    return x * x


print(opti_tools.evalutate_at_10(square))

print(opti_tools.solve_euler(0.0, 10.0, 10000, 1.0, v))
