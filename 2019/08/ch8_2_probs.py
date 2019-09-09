import numpy as np


def improved_euler(f, y0, h, ts, method='improved'):

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

        if method == 'improved':

            # Compute the values needed in the improved Euler scheme.
            k1 = f(t, yn)
            k2 = f(t + h, yn + h * k1)
            ynp1 = yn + h / 2.0 * (k1 + k2)

        elif method == 'modified':

            k1 = f(t + 0.5 * h, yn + 0.5 * h * f(t, yn))
            ynp1 = yn + h * k1

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
    h  = 0.0125
    ts = [0.1, 0.2, 0.3, 0.4]

    # Our f(t, y) for this problem.
    def f(t, y):
        return 3 + t - y

    improved_euler(f, y0, h, ts)

def prob10():

    # Starting values.
    y0 = 1
    h  = 0.0125
    ts = [0.5, 1.0, 1.5, 2.0]

    def f(t, y):
        return 2 * t + np.exp(-t * y)

    improved_euler(f, y0, h, ts)

def prob23():

    # Starting values.
    y0 = 1
    h  = 0.05
    ts = [0.1, 0.2, 0.3, 0.4]

    # Our f(t, y) for this problem. Same as problem 1 and results are to be compared.
    def f(t, y):
        return 3 + t - y

    improved_euler(f, y0, h, ts, 'modified')
