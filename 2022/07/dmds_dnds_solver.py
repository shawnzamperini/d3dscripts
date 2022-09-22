# Use sympy to solve for equations of dM/ds and dn/ds.
import sympy as sm
from sympy import symbols


# Define symbols.
dMds, M, S, n, cs, A, Te, m, dAds, dnds = symbols("dMds, M, S, n, cs, A, Te, m, dAds, dnds", real=True)

# Declare RHS's as expressions.
rhs1 = S/(n*cs) - dnds*M/n - M/A*dAds
rhs2 = -m*n*M*cs**2*dMds/(2*Te) - m*M*cs*S/(2*Te)

# Solver looks for equations of the form, e.g., dM/ds - RHS = 0.
result = sm.solve([dMds - rhs1, dnds - rhs2], [dMds, dnds])

# Print out, easy to copy/paste.
for key in result.keys():
    print("{} = {}".format(key, result[key]))
