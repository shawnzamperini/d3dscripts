import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle
from scipy.interpolate import interp2d
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from scipy.signal import savgol_filter


# Match each RCP file with the EFIT data.
efit_file = {"MP187103_1":"187103_1620", "MP187103_2":"187103_3410",
             "MP187104_1":"187104_1620", "MP187104_2":"187104_3410",
             "MP187105_1":"187105_1620", "MP187105_2":"187105_3410",
             "MP187106_1":"187106_1620", "MP187106_2":"187106_3410",
             "MP187107_1":"187107_1620", "MP187107_2":"187107_3410",
             "MP187108_1":"187108_1620", "MP187108_2":"187108_3410",
             "MP187109_1":"187109_1620", "MP187109_2":"187109_3410",
             "MP187110_1":"187110_1620", "MP187110_2":"187110_3410",
             "MP187111_1":"187111_1620", "MP187111_2":"187111_3410"}

pinjs = {187103:2500, 187104:2500, 187105:2500, 187106:4400, 187107:4400,
        187108:4400, 187109:3300, 187110:4400, 187111:3300}

rcp_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/rcp_data/rcp_master.xlsx"
#rcp_dfs = pd.read_excel(rcp_path, sheet_name=None)
rcp_xl = pd.ExcelFile(rcp_path)
rcp_names = rcp_xl.sheet_names

def exp_fit(x, a, b, c):
    return a * np.exp(-x * b) + c

def load_pickle(path):
    with open(path, "rb") as f:
        var = pickle.load(f)
    return var

plot_data = {}
for rcp_name in rcp_names:

    # Do things a little differently so we can get the sheet names.
    rcp_df = rcp_xl.parse(rcp_name)

    # Not all data has floating potential yet, so head those off here with a
    # continue statement.
    try:
        rcp_df["Vf1(V)"]
    except:
        continue

    # Load in the respective EFIT data.
    root = "/mnt/c/Users/Shawn/Documents/d3d_work/gfile_data_for_rcp/"
    try:
        r_efit  = load_pickle(root + efit_file[rcp_name] + "_r")
        z_efit  = load_pickle(root + efit_file[rcp_name] + "_z")
        bp_efit = load_pickle(root + efit_file[rcp_name] + "_bp")
        bt_efit = load_pickle(root + efit_file[rcp_name] + "_bt")
    except:
        print("{}: Error loading pickle data.".format(rcp_name))
        continue

    # Interpolations along the (R, Z) values of the RCP to determine Bt and Bp
    # at each measurement location.
    #R_efit, Z_efit = np.meshgrid(r_efit, z_efit)
    #f_bp = Rbf(R_efit, Z_efit, bp_efit)
    f_bp = interp2d(r_efit, z_efit, bp_efit)
    f_bt = interp2d(r_efit, z_efit, bt_efit)

    # Remove all data beyond R = 2.35m (2.36m is where I think it starts going
    # outside the vessel).
    #print("Removing data beyond 2.35m...")
    #mask = rcp_df["R(cm)"] / 100 < 2.35
    #rcp_df = rcp_df[mask]

    rho = rcp_df["Rho"]
    if rcp_name[:2] == "MP":
        r = rcp_df["R(cm)"] / 100
        z = np.full(len(rcp_df), -0.188)
    elif rcp_name[:2] == "XP":
        r = np.full(len(rcp_df), 1.493)
        z = rcp_df["Z(cm)"] / 100
    else:
        print("Error: Did not identify probe name {}".format(rcp_name))

    bp = f_bp(r, z)[0]
    bt = f_bt(r, z)[0]
    b  = np.sqrt(bp**2 + bt**2)

    # Average potential is that of these two.
    rcp_df["Vf_avg"] = (rcp_df["Vf1(V)"] + rcp_df["Vf2(V)"]) / 2.0 + 3 * rcp_df["Te(eV)"]
    vf_exp = rcp_df["Vf_avg"]

    # Just a basic third order fit should be good enough.
    poly_order = 3
    poly3_vf = np.polynomial.Polynomial(poly_order).fit(r, rcp_df["Vf_avg"], poly_order)
    vf_fit = poly3_vf(r)
    er_fit = -poly3_vf.deriv()(r)

    #sg_window = 7; sg_poly = 2
    #vf_fit = savgol_filter(vf_exp, sg_window, sg_poly)
    #er_fit = -np.gradient(vf_fit, r)

    # Ultimately we are interested in the ion pressure (temperature) gradient,
    # not the electron.
    te_mult = 1.0
    print("Multiplying Te by {:.2f}".format(te_mult))

    # Nothing stopping us from calculating the rafdial pressure gradient either.
    # Prolly best to do exp. fits maybe to ne, Te first though?
    te_exp = rcp_df["Te(eV)"] * te_mult
    #popt_te, pcov_te = curve_fit(exp_fit, r, te_exp, maxfev=5000, p0=(0.01, 0.1, 1.0))
    #te_fit = exp_fit(r, *popt_te)
    ne_exp = rcp_df["Ne (E18 m-3)"]
    #popt_ne, pcov_ne = curve_fit(exp_fit, r, ne_exp, maxfev=5000, p0=(0.01, 0.1, 1.0))
    #ne_fit = exp_fit(r, *popt_ne)
    #pe_exp = ne_exp * 1e18 * te_exp
    #pe_fit = ne_fit * 1e18 * te_fit

    # Or maybe a polynomial again...
    poly3_te = np.polynomial.Polynomial(poly_order).fit(r, te_exp, poly_order)
    te_fit = poly3_te(r)
    poly3_ne = np.polynomial.Polynomial(poly_order).fit(r, ne_exp, poly_order)
    ne_fit = poly3_ne(r)

    #te_fit = savgol_filter(te_exp, sg_window, sg_poly)
    #ne_fit = savgol_filter(ne_exp, sg_window, sg_poly)

    e = 1.609e-19
    pe_exp = ne_exp * 1e18 * te_exp * e
    pe_fit = ne_fit * 1e18 * te_fit * e
    pe_fit_grad = (poly3_te(r) * poly3_ne.deriv()(r) + poly3_ne(r) * poly3_te.deriv()(r)) * 1e18 * e
    #pe_fit_grad = (te_fit * np.gradient(ne_fit, r)) + (ne_fit * np.gradient(te_fit, r)) * 1e18 * e

    # Outer midplane approximation. Need to update with the real one that
    # accounts for angle eventually.
    a = 0.593
    v_comb = er_fit / bp - 2 * (a / r) * (b / bp) * pe_fit_grad / (e * ne_fit * 1e18 * b)
    v_flow = rcp_df["Vflow (km/s)"] * 1000

    # PF flow midplane approx.
    v_pf = 2 * (a / r) * (b / bp) * (er_fit / b - pe_fit_grad / (e * ne_fit * 1e18 * b))

    shot = int(rcp_name[2:8])
    pinj = np.array(list(pinjs.values()))
    color_code = (pinjs[shot] - pinj.min()) / (pinj.max() - pinj.min())

    plot_data[rcp_name] = {"R": r, "Vf_exp":vf_exp,
      "Vf_fit":vf_fit, "Er_fit":er_fit, "Mach":rcp_df["Machn"],
      "Te_exp":te_exp, "Te_fit":te_fit, "ne_exp":ne_exp, "ne_fit":ne_fit,
      "pe_exp":pe_exp, "pe_fit":pe_fit, "pe_fit_grad":pe_fit_grad,
      "v_comb":v_comb, "v_flow":v_flow, "shot":shot, "color_code":color_code,
      "rho":rho, "v_pf":v_pf, "bp":bp, "bt":bt, "b":b}

# Get min and max rho values.
rho_min = 99999; rho_max = 0
for rcp_name, data in plot_data.items():
    if data["rho"].min() < rho_min:
        rho_min = data["rho"].min()
    if data["rho"].max() > rho_max:
        rho_max = data["rho"].max()

ml_df = pd.DataFrame()

# 9x9 grids of plots.
count = 0
plt.close("all")
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
for i in range(0, len(plot_data.keys())):
    if i % 9 == 0:

        # Show previous figure if we aren't at the start of the loop.
        if i > 0:
            fig.tight_layout()
            fig.show()

        fig, axs = plt.subplots(3, 3, figsize=(10, 8))
        axs = axs.flatten()
        count = 0

    # Pull out some data.
    keys = list(plot_data.keys())
    r      = plot_data[keys[i]]["R"].values
    rho    = plot_data[keys[i]]["rho"].values
    vf_exp = plot_data[keys[i]]["Vf_exp"].values
    vf_fit = plot_data[keys[i]]["Vf_fit"]
    er_fit = plot_data[keys[i]]["Er_fit"]
    mach   = plot_data[keys[i]]["Mach"].values
    te_exp = plot_data[keys[i]]["Te_exp"].values
    te_fit = plot_data[keys[i]]["Te_fit"]
    ne_exp = plot_data[keys[i]]["ne_exp"].values
    ne_fit = plot_data[keys[i]]["ne_fit"]
    pe_exp = plot_data[keys[i]]["pe_exp"].values
    pe_fit = plot_data[keys[i]]["pe_fit"]
    pe_fit_grad = plot_data[keys[i]]["pe_fit_grad"]
    v_comb = plot_data[keys[i]]["v_comb"]
    v_flow = plot_data[keys[i]]["v_flow"].values
    v_pf   = plot_data[keys[i]]["v_pf"]

    # Organize into our ML DataFrame.
    tmp = dict(zip(["r", "rho", "vf_fit", "er_fit", "mach", "te_fit", "ne_fit",
      "pe_fit", "pe_fit_grad", "v_comb", "v_flow"], [r, rho, vf_fit, er_fit,
      mach, te_fit, ne_fit, pe_fit, pe_fit_grad, v_comb, v_flow]))
    df = pd.DataFrame(tmp)
    ml_df = ml_df.append(df)

    # Can change plot choices here.
    y_exp = vf_exp; y_fit = vf_fit; y_sec = er_fit

    axs[count].scatter(r, y_exp, c="tab:red", edgecolors="k")
    axs[count].plot(r, y_fit, color="tab:red")
    axs[count].text(0.05, 0.80, keys[i], transform=axs[count].transAxes)
    ax_sec = axs[count].twinx()
    ax_sec.plot(r, y_sec, color="tab:purple")

    count += 1

    if i % 2 == 999:
        pass
    else:
        ax2.scatter(er_fit[3:-3], mach[3:-3], edgecolors="k", label=keys[i])

        cmap = plt.get_cmap('Reds')
        #color_code = (rho - rho_min) / (rho_max - rho_min)
        shot = int(keys[i][2:8])
        color_code = pinjs[shot]
        pinj_max = max(list(pinjs.values())); pinj_min = min(list(pinjs.values()))
        color_code = (pinjs[shot] - pinj_min) / (pinj_max - pinj_min)
        ax3.scatter(v_comb, v_flow, edgecolors="k", label=keys[i], color=cmap(color_code))

        ax4.scatter(te_fit, v_flow - v_comb, edgecolors="k", label=keys[i])

fig.tight_layout()
fig.show()

ax2.legend(ncol=2)
ax2.grid()
ax2.set_xlabel("Er")
ax2.set_ylabel("Mach")
fig2.tight_layout()
fig2.show()

ax3.plot(ax3.get_ylim(), ax3.get_ylim(), "k--")
ax3.set_xlim(ax3.get_ylim())
ax3.legend(ncol=2)
ax3.grid()
ax3.set_xlabel("V Comb")
ax3.set_ylabel("Measured Flow")
fig3.tight_layout()
fig3.show()

ax4.legend(ncol=2)
ax4.set_xlabel("Te")
ax4.set_ylabel("v_flow - v_comb")
fig4.tight_layout()
fig4.show()

# Some ML afterwards to see if we can pick up on what's causing the most discrepancy.
X = ml_df[["r", "ne_fit", "te_fit", "vf_fit", "er_fit"]]
#X = ml_df[["r", "v_comb"]]
y = ml_df["v_flow"]
X_train, X_test, y_train, y_test = train_test_split(X, y)
scaler = StandardScaler().fit(X_train)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)
model = RandomForestRegressor(n_estimators=1000).fit(X_train_scaled, y_train)
#model = LinearRegression().fit(X_train_scaled, y_train)
print("Train: {:.3f}".format(model.score(X_train_scaled, y_train)))
print("Test:  {:.3f}".format(model.score(X_test_scaled, y_test)))

fig5, ax5 = plt.subplots()
ax5.scatter(model.predict(X_test_scaled), y_test, edgecolors="k")
ax5.set_xlabel("Model")
ax5.set_ylabel("v_flow")
ax5.plot(ax5.get_ylim(), ax5.get_ylim(), "k--")
ax5.set_xlim(ax5.get_ylim())
fig5.tight_layout()
fig5.show()
