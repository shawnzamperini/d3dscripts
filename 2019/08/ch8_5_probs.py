import numpy as np


def solve_system(fs, xy0, h, ts, method):

    # Start at t=0.
    t = 0; n = 0

    # z will represent the array (vector) of the solution, where the two
    # components are x and y.
    #
    # z = (x)
    #     (y)
    zn = xy0

    # Print the header.
    print("{:<5} {:<5} {:<7} {:<7}".format("n", "t", "xn", "yn"))

    # These are arrays to hold values for the two methods. Euler only uses k1.
    k1 = np.zeros((2)); k2 = np.zeros((2)); k3 = np.zeros((2)); k4 = np.zeros((2));

    # Assuming ts are in order, go up until you hit the highest t.
    while t <= ts[-1]:

        # Round the t off to prevent floating point errors.
        t    = np.round(t, 5)
        tnp1 = np.round(t + h, 5)

        if method == 'euler':

            # Evaluate the function for each diff. eq.
            k1[0] = fs[0](t, zn[0], zn[1])
            k1[1] = fs[1](t, zn[0], zn[1])

            # Euler method applied just doing array style math.
            znp1 = zn + h * k1

        elif method == 'runge_kutta':

            # Solve for the values to be put into the RK scheme. Would be better
            # to make the functions accept zn as an array.
            k1[0] = fs[0](t, zn[0], zn[1])
            k1[1] = fs[1](t, zn[0], zn[1])

            k2[0] = fs[0](t + 0.5 * h, zn[0] + 0.5 * h * k1[0], zn[1] + 0.5 * h * k1[1])
            k2[1] = fs[1](t + 0.5 * h, zn[0] + 0.5 * h * k1[0], zn[1] + 0.5 * h * k1[1])

            k3[0] = fs[0](t + 0.5 * h, zn[0] + 0.5 * h * k2[0], zn[1] + 0.5 * h * k2[1])
            k3[1] = fs[1](t + 0.5 * h, zn[0] + 0.5 * h * k2[0], zn[1] + 0.5 * h * k2[1])

            k4[0] = fs[0](t + h, zn[0] + h * k3[0], zn[1] + h * k3[1])
            k4[1] = fs[1](t + h, zn[0] + h * k3[0], zn[1] + h * k3[1])

            # Solve for z_n+1.
            znp1 = zn + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        # If we are at one of the ts where we want to print the results.
        if t in ts:
            print("{:<5} {:<5} {:<7.6} {:<7.6}".format(n, t, zn[0], zn[1]))

        # Update value for next iterations.
        n = n + 1
        t = tnp1
        zn = znp1

def prob1():

    # Inital values and inputs.
    x0 = 1
    y0 = 0
    h  = 0.2
    ts = np.array([0.2, 0.4, 0.6, 0.8, 1.0])

    # The equation for x'.
    def f1(t, x, y):
        return x + y + t

    # The equation for y'.
    def f2(t, x, y):
        return 4 * x - 2 * y

    # Pass the functions as an array (i.e. vector), as well as the initial conditions.
    fs  = np.array([f1, f2])
    xy0 = np.array([x0, y0])

    solve_system(fs, xy0, h, ts, 'runge_kutta')
