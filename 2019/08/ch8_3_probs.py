import numpy as np


def runge_kutta(f, y0, h, ts):

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

        k1 = f(t,           yn)
        k2 = f(t + 0.5 * h, yn + 0.5 * h * k1)
        k3 = f(t + 0.5 * h, yn + 0.5 * h * k2)
        k4 = f(t + h,       yn + h * k3)
        ynp1 = yn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

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

    # Starting values.
    y0 = 1
    h  = 0.1
    ts = [0.1, 0.2, 0.3, 0.4]

    # Our f(t, y) for this problem.
    def f(t, y):
        return 3 + t - y

    runge_kutta(f, y0, h, ts)

def prob5():

    # Starting values.
    y0 = 0.5
    h  = 0.1
    ts = [0.1, 0.2, 0.3, 0.4]

    # Our f(t, y) for this problem.
    def f(t, y):
        return (y**2 + 2 * t * y) / (3 + t**2)

    runge_kutta(f, y0, h, ts)
