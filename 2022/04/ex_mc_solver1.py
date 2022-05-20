# This script is a super basic demonstration of how to "solve" the equation:
#       dx/dt = sin(t) + t^3            x(0) = -1
# via a Monte-Carlo typie method. The solution is
#      x(t) = -cos(t) + t^4 / 4
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
from scipy.interpolate import interp1d


N = 500  # Number of points in each solution.
domain = [0, 10]
NN = 10  # Number of solutions to generate.
smooth = None  # One of savgol or poly or None
smooth_poly = 4
width_mult = 0.5
NNN = 100  # Number of iterations, where each iteration uses the best previous solution.
iterate_width = True
grade_type = "residuals"  # One of residuals, rmse

def fps(x, y, i):
    """
    Implement a five-point stencil for smoothing the derivative some.
    https://en.wikipedia.org/wiki/Five-point_stencil

    Assumes evenly spaced data.
    f0 = f(x)  -  Not used
    f1 = f(x+h)
    f2 = f(x+2h)
    f3 = f(x-h)
    f4 = f(x-2h)

    --- Input ---
    x: The x data
    y: The y data
    i: The index to get the derivative at
    """

    # Special case for the edges.
    if i < 2:
        return (y[i+1] - y[i]) / (x[i+1] - x[i])
    elif i > len(x) - 3:
        return (y[i] - y[i-1]) / (x[i] - x[i-1])

    else:
        h = x[i] - x[i-1]
        f0 = y[i]
        f1 = y[i+1]
        f2 = y[i+2]
        f3 = y[i-1]
        f4 = y[i-2]
        return (-f2 + 8*f1 - 8*f3 + f4) / (12 * h)

# Choose t randomly.
#t = domain[0] + (domain[1] - domain[0]) * random.random()
t_vals = np.linspace(domain[0], domain[1], N)

# Using the ODE, fill out values of the dependent variable.
dxdt_vals = np.sin(t_vals) + np.power(t_vals, 3)

best_solutions = []
prev_best_grade = 1e10
best_t = np.zeros(N)
best_x = np.zeros(N)
best_res = np.full(N, 1e10)

for k in tqdm(range(0, NNN)):
    all_results = []
    for j in range(0, NN):
        results = {"t":np.zeros(N), "x":np.zeros(N), "dxdt":np.zeros(N)}
        #for i in range(0, N):

        # Store results.
        results["t"] = t_vals
        results["dxdt"] = dxdt_vals

        # Sort the values according to t.
        #sort_idx = np.argsort(results["t"])
        #results["t"] = results["t"][sort_idx]
        #results["dxdt"] = results["dxdt"][sort_idx]

        # Construct the "guess" for x. A reasonable guess would be t^4 / 4 here (in
        # fact a very good "guess"). To add variation we will sample from it via a
        # normal distribution.
        if k == 0:
            #guess_x = np.power(results["t"], 4) / 4
            results["guess_x"] = np.power(results["t"], 3) / 3
        else:
            results["guess_x"] = f_best(results["t"])

        # Now, at each t, use the two bounding points to get the point-slope equation
        # of dxdt = mt + b. x at this t location is then x = 1/2 * mt^2 + bt + c
        # with the boundary condition x(0) = -1 to determine c.
        slopes = []
        for i in range(0, N):

            # The first and last points we just have to use two points.
            if i == 0:
                ts = results["t"][0:2]
                dxdts = results["dxdt"][0:2]
            elif i == N-1:
                ts = results["t"][-2:]
                dxdts = results["dxdt"][-2:]
            else:
                ts = results["t"][i-1:i+2]
                dxdts = results["dxdt"][i-1:i+2]

            t = results["t"][i]

            # Do linear fit, but then what's the BC for this "temporary" integral of
            # the linear fit to determine c? Let's try with an intial guess on what
            # x could be.
            #m, b = np.polyfit(ts, dxdts, 1)
            center = results["guess_x"][i]
            if k == 0 or not iterate_width:
                width = np.abs(center * width_mult)
            else:
                if grade_type == "residuals":
                    width = widths[i]
                else:
                    width = prev_best_grade
            x_ran = np.random.normal(center, width)
            #c = x_ran - 0.5 * m * t**2 - b * t
            #x = 0.5 * m * t**2 + b * t + c
            x = x_ran
            results["x"][i] = x
            #slopes.append([ts, m*ts+b])

        # Offer a smoothed result.
        if smooth == "savgol":
            window = N / 10
            if window % 2 == 0:
                window += 1
            window = int(window)
            results["x_smooth"] = savgol_filter(results["x"], window, 2)
        elif smooth == "poly":
            z = np.polyfit(results["t"], results["x"], smooth_poly)
            p = np.poly1d(z)
            results["x_smooth"] = p(results["t"])
        else:
            results["x_smooth"] = results["x"]

        all_results.append(results)


    # Now for every result, go through and grade the solutions by
    grades = []
    for i in range(0, len(all_results)):
        results = all_results[i]
        t = results["t"]

        # Need to recalculate dxdt from the solved x-values.
        #dxdt = np.gradient(results["x_smooth"], t)
        dxdt = []
        for m in range(0, len(t)):
            dxdt.append(fps(t, results["x_smooth"], m))
        lhs = np.array(dxdt)
        rhs = np.sin(t) + np.power(t, 3)

        if grade_type == "rmse":
            # Is using a savgol we don't want to grade the solution that gets
            # extrapolated outside the window.
            if smooth == "savgol":
                lhs = lhs[window:-window]
                rhs = rhs[window:-window]
            grades.append(np.sqrt(mean_squared_error(lhs, rhs)))
        elif grade_type == "residuals":
            grades.append(np.square(lhs-rhs))


    # Identify the best solution. This solution will feed
    # in as the guess for the next iteration only if it's better than the previous
    # best solution.
    if grade_type == "rmse":
        best = np.argmin(grades)
        print("Best RMSE: {:.3f}".format(grades[best]))
        if grades[best] < prev_best_grade:
            f_best = interp1d(all_results[best]["t"], all_results[best]["x_smooth"], fill_value="extrapolate")
            best_solutions.append(all_results[best])
            prev_best_grade = grades[best]

    # In this weird case just grab the data point from all the curves at each
    # location that has the lowest residual.
    elif grade_type == "residuals":
        for i in range(0, NN):
            result = all_results[i]
            for ii in range(0, N):
                res = grades[i][ii]
                if res < best_res[ii]:
                    best_t[ii] = result["t"][ii]
                    best_x[ii] = result["x_smooth"][ii]
                    best_res[ii] = res

        best_solutions.append({"t":np.copy(best_t), "x_smooth":np.copy(best_x)})
        f_best = interp1d(best_t, best_x, fill_value="extrapolate")
        best = -1

        # Likely good idea to make the widths the lingering residuals.
        widths = best_res


ans_t = np.linspace(domain[0], domain[1], 100)
ans_x = -np.cos(ans_t) + np.power(ans_t, 4) / 4
fig, ax = plt.subplots()
alphas = np.linspace(0.1, 1, len(best_solutions))
if NNN < 100:
    for i in range(0, len(best_solutions)):
        t = best_solutions[i]["t"]
        x = best_solutions[i]["x_smooth"]
        ax.plot(t, x, alpha=alphas[i], color="r")
        #ax.scatter(results["t"], results["x"])
ax.plot(ans_t, ans_x, lw=3, color="k")
if grade_type == "residuals":
    ax.plot(best_t, best_x, lw=3, color="r")
else:
    ax.plot(all_results[best]["t"], all_results[best]["x_smooth"], lw=3, color="r")
ax.set_xlabel("t")
ax.set_ylabel("x")
fig.tight_layout()
fig.show()
