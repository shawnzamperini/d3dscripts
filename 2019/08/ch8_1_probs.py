import numpy as np
import mpmath as mp
from sympy.solvers import solve
from sympy import Symbol



def euler(method, f, y0, h, ts, solved_ynp1=None):

    # Start at t=0.
    t = 0; n = 0

    # Start at the intial value (y(t=0) in all these examples).
    yn = y0

    # Print the header.
    print("{:<5} {:<5} {:<5}".format("n", "t", "yn"))

    # Assuming ts are in order, go up until you hit the highest t.
    while t <= ts[-1]:

        # Round the t off to prevent floating point errors.
        t    = np.round(t, 5)
        tnp1 = np.round(t + h, 5)

        if method == 'forward':

            # The next estimate, y_n+1, is given by y_n+1 = yn + hf_n.
            ynp1 = yn + h * f(t, yn)

        elif method == 'reverse':

            # Need to have solution of the implicit equation for ynp1.
            if solved_ynp1 is None:
                print("Error: Must supply a solved_ynp1 for Reverse Euler.")
                break
            else:
                #print(yn)
                ynp1 = solved_ynp1(yn, h, tnp1)

        else:
            print("Error: Incorrect method entry.")
            break

        # If we are at one of the ts where we want to print the results.
        if t in ts:
            print("{:<5} {:<5} {:<5}".format(n, t, yn))

        # Assign t to t_n+1 for next iteration. Increase index. Assign yn to
        # ynp1 for next round.
        t  = tnp1
        n  = n + 1
        yn = ynp1

    # Put a blank line at the end.
    print("")


def prob1():
    ts = [0.1, 0.2, 0.3, 0.4]

    def f(t, y):
        return 3 + t - y

    def solved_ynp1(yn, h, tnp1):
        return 1/(1+h) * (yn + 3*h + h*tnp1)

    print("Forward Euler, h = 0.05")
    euler("forward", f=f, y0=1, h=0.05,  ts=ts)
    print("Forward Euler, h = 0.025")
    euler("forward", f=f, y0=1, h=0.025, ts=ts)
    print("Reverse Euler, h = 0.05")
    euler("reverse", f=f, y0=1, h=0.05,  ts=ts, solved_ynp1=solved_ynp1)
    print("Reverse Euler, h = 0.025")
    euler("reverse", f=f, y0=1, h=0.025, ts=ts, solved_ynp1=solved_ynp1)

def prob6():
    ts = [0.1, 0.2, 0.3, 0.4]

    def f(t, y):
        return (t**2 - y**2) * np.sin(y)

    def solved_ynp1(yn, h, tnp1):
        return None

    print("Forward Euler, h = 0.05")
    euler("forward", f=f, y0=-1, h=0.05,  ts=ts)
    print("Forward Euler, h = 0.025")
    euler("forward", f=f, y0=-1, h=0.025, ts=ts)
    print("Reverse Euler, h = 0.05")
    euler("reverse", f=f, y0=-1, h=0.05,  ts=ts, solved_ynp1=solved_ynp1)
    print("Reverse Euler, h = 0.025")
    euler("reverse", f=f, y0=-1, h=0.025, ts=ts, solved_ynp1=solved_ynp1)

def prob11():
    ts = [0.5, 1.0, 1.5, 2.0]

    def f(t, y):
        return (4-t*y)/(1+y**2)

    def solved_ynp1(yn, h, tnp1):
        y = Symbol('y')
        solve(yn + h * (4 - tnp1 * y) / (1 + y**2), y)

    print("Forward Euler, h = 0.05")
    euler("forward", f=f, y0=-2, h=0.025,  ts=ts)
    print("Forward Euler, h = 0.025")
    euler("forward", f=f, y0=-2, h=0.0125, ts=ts)
    print("Reverse Euler, h = 0.05")
    euler("reverse", f=f, y0=-2, h=0.025,  ts=ts, solved_ynp1=solved_ynp1)
    print("Reverse Euler, h = 0.025")
    euler("reverse", f=f, y0=-2, h=0.0125, ts=ts, solved_ynp1=solved_ynp1)
