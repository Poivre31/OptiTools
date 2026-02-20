"""Mon module de controle optimal"""

from collections.abc import Callable


def solve_euler(t0: float, tf: float, N: int, x0: float, v: Callable[[float, float], float]) -> float:
    """Calcule x(tf) pour l'état initial x0 sous la contrainte dx/dt=v(t,x)"""

class solver:
    def __init__(self, arg: int, /) -> None: ...

    def solve(self, arg: float, /) -> float: ...
