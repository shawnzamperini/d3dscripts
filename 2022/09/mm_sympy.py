from sympy.solvers import solve
from sympy import Symbol, solveset
from sympy.solvers.solveset import linsolve


x1 = Symbol("x1")
x2 = Symbol("x2")
x3 = Symbol("x3")
x4 = Symbol("x4")
x5 = Symbol("x5")
x6 = Symbol("x6")
x7 = Symbol("x7")
x8 = Symbol("x8")

A1 = Symbol("A1")
A2 = Symbol("A2")
A3 = Symbol("A3")
A4 = Symbol("A4")
A5 = Symbol("A5")
A6 = Symbol("A6")
A7 = Symbol("A7")
A8 = Symbol("A8")
A9 = Symbol("A9")
A10 = Symbol("A10")
A11 = Symbol("A11")
A12 = Symbol("A12")
A13 = Symbol("A13")
A14 = Symbol("A14")
A15 = Symbol("A15")
A16 = Symbol("A16")

R = Symbol("R")

eqs = [
    A1 + A2 + x6 * A3 + x5 * A13 - x1,
    A4 + A5 + x6 * A6 + x5 * A14 - x2,
    A7 + A8 * x6 * A9 + x5 * A15 - x3,
    A10 + A11 * x6 * A12 + x5 * A16 - x4,
    (1 - x7 - x8) * x1 + x8 * x3 - x5,
    (1 - x7 - x8) * x2 + x7 * x4 - x6,
    (1 - R) * x6 / x4 - x7,
    (1 - (1 - R) * x6 / x4) * (x2 - x1) / (x3 + x2 - x1) - x8]

solve(eqs, (x1, x2, x3, x4, x5, x6, x7, x8))
