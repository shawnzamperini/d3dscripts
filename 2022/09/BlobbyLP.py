# This class represents a workflow to use a Monte Carlo algorithm to
# interpretively model LP target ne profiles.
import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import pandas as pd


# Need this path for ThomsonClass.
sys.path.append("/Users/zamperini/github/utk-fusion/tools/")

class BlobbyLP:

    def __init__(self):
        pass

    def calc_te_params(self, shot, tmin, tmax, tmany):
        """
        Pull Thomson Scattering data and provide an interface for fitting the
        Te data so that we can obtain Te,sep and lambda_Te.

        Inputs
        shot (int): Shot number.
        tmin (float): Minimum time for averaging.
        tmax (float): Maximum time for averaging.
        tmany (int): How many time slices (gfiles) to load for averaging.

        """

        from ThomsonClass import ThomsonClass
        from scipy.optimize import curve_fit

        # Load the TS data.
        ts = ThomsonClass(shot, 'core')
        ts.load_ts(tunnel=False)
        ts.map_to_efit(np.linspace(tmin, tmax, tmany), tree="EFIT01")

        # Pull out the arrays.
        r  = ts.avg_omp['RminRsep_omp'] * 100
        te = ts.avg_omp['Te_omp']
        ne = ts.avg_omp['ne_omp'] * 10**(-18)
        r_err  = ts.avg_omp['RminRsep_omp_err'] * 100  # m to cm
        te_err = ts.avg_omp['Te_omp_err']
        ne_err = ts.avg_omp['ne_omp_err'] * 10**(-18) # m-3 to 10^-18 m-3

        # Plot with chord numbers next to each so you can point out the bad ones.
        plt.ion()
        fig = plt.figure(figsize=(10,5))
        ax_te = fig.add_subplot(121)
        ax_ne = fig.add_subplot(122)
        ax_te.errorbar(r, te, te_err, r_err, fmt='k.', ms=8)
        ax_ne.errorbar(r, ne, ne_err, r_err, fmt='k.', ms=8)
        ax_te.set_xlabel('R-Rsep OMP(cm)', fontsize=14)
        ax_ne.set_xlabel('R-Rsep OMP(cm)', fontsize=14)
        ax_te.set_ylabel('Te (eV)', fontsize=14)
        ax_ne.set_ylabel(r'$\mathdefault{ne\ (10^{18}\ m^{-3})}$', fontsize=14)
        for i, chord in enumerate(np.arange(0, len(r))):
            ax_te.annotate(chord, (r[i], te[i]))
            ax_ne.annotate(chord, (r[i], ne[i]))
        fig.tight_layout()
        fig.show()

        while True:

            # Choose chords for fitting. Just start with the first three inside the sep
            # and the rest outwards.
            num_in_include = 3
            first_out = np.where(r > 0)[0][-1]
            include = num_in_include + first_out

            # Zoom into the fitting region on the plots.
            ax_te.set_xlim(1.1 * r[:include].min(), 1.1 * r[:include].max())
            ax_te.set_ylim(1.1 * te[:include].min(), 1.1 * te[:include].max())
            ax_ne.set_xlim(1.1 * r[:include].min(), 1.1 * r[:include].max())
            ax_ne.set_ylim(1.1 * ne[:include].min(), 1.1 * ne[:include].max())

            # Ask if any extra chords should be excluded.
            print('\nChords for fitting: ', end=''); print(*np.arange(0, len(r))[:include])
            exclude = input('Chords to exclude (separated by commas, press enter if none): ').split(',')
            exclude = np.array(exclude, dtype=int)

            if exclude != ['']:
                r_tofit  = np.delete(r[:include],  exclude)
                te_tofit = np.delete(te[:include], exclude)
                ne_tofit = np.delete(ne[:include], exclude)
            else:
                r_tofit  = np.array(r[:include])
                te_tofit = np.array(te[:include])
                ne_tofit = np.array(ne[:include])

            def exp_fit(x, a, b):
                return a * np.exp(-x*b)

            # Perform the fitting.
            popt_te, pcov_te = curve_fit(exp_fit, r_tofit, te_tofit, maxfev=10000)
            popt_ne, pcov_ne = curve_fit(exp_fit, r_tofit, ne_tofit, maxfev=10000)

            r_fit = np.linspace(r_tofit.min(), r_tofit.max(), 50)
            te_fit = exp_fit(r_fit, *popt_te)
            ne_fit = exp_fit(r_fit, *popt_ne)

            # Plot the fit.
            ax_te.plot(r_fit, te_fit, 'k--', lw=5)
            ax_ne.plot(r_fit, ne_fit, 'k--', lw=5)
            #fig.show()

            print('Te: {:.2f} * exp(-r / {:.2f})'.format(popt_te[0], 1/popt_te[1]))
            print('ne: {:.2f} * exp(-r / {:.2f})'.format(popt_ne[0], 1/popt_ne[1]))

            print('\nTe falloff (cm): {:.2f}'.format(1 / popt_te[1]))
            print('ne falloff (cm): {:.2f}'.format(1 / popt_ne[1]))

            ans = input('Try fitting again? (y/n): ')
            if ans == 'n':
                break

        # Return the fit values.
        tesep = popt_te[0]
        nesep = popt_ne[0]
        lambda_te = 1 / popt_te[1]
        lambda_ne = 1 / popt_ne[1]
        return {"tesep":tesep, "nesep":nesep, "lambda_te":lambda_te,
            "lambda_ne":lambda_ne}

    def calc_flux_expansion(self, shot, time, gfile_pickle_path):
        """
        Calculate the flux expansion along the target using the EFIT. This
        returns a dictionary of the flux expansion FROM THE X-POINT in relation
        to the target since we assume blobs only travel the X-point distance
        to get to the target. A flux expansion of 1.0 is at the target, while
        greater than 1 is towards/at the X-point.

        Inputs
        shot (int): Shot number.
        gfile_pickle_path (str): Path the the pickled gfile.
        """

        from scipy.interpolate import RectBivariateSpline, interp1d, Rbf

        # Load gfile.
        with open(gfile_pickle_path, "rb") as f:
            gfile = pickle.load(f)

        # Need the location of the strike point (RVSOUT).
        from gadata import gadata
        import MDSplus
        conn = MDSplus.Connection("atlas.gat.com")
        gaobj = gadata("RVSOUT", shot, connection=conn)
        times = gaobj.xdata
        rvsouts = gaobj.zdata
        idx = np.argmin(np.abs(times-time))
        rvsout = rvsouts[idx]

        # Get the total and poloidal magnetic fields along the shelf. Can save
        # time by restricting the region to just the shelf region.
        Rs, Zs = np.meshgrid(gfile["R"], gfile["Z"])
        #f_bt = RectBivariateSpline(gfile["R"], gfile["Z"], gfile["Bt"])
        #f_bp = RectBivariateSpline(gfile["R"], gfile["Z"], gfile["Bp"])

        # Get coordinates of a line for the shelf (easy) and OMP (harder).
        shelf_r1 = 1.372
        shelf_r2 = 1.678
        shelf_z  = -1.25

        # REPLACING with the X-point, so that we are calculating the flux
        # expansion between the OMP and the X-point.
        #gaobj = gadata("RXPT1", shot, connection=conn)
        #times = gaobj.xdata
        #rxpts = gaobj.zdata
        #idx = np.argmin(np.abs(times-time))
        #shelf_r1 = rxpts[idx]
        #shelf_r2 = 1.85
        #gaobj = gadata("ZXPT1", shot, connection=conn)
        #times = gaobj.xdata
        #zxpts = gaobj.zdata
        #idx = np.argmin(np.abs(times-time))
        #shelf_z = zxpts[idx]


        # OMP. First get an interpolation of the separatrix, just the outboard
        # side. Then see what the R value of the sepratrix is at the OMP. Then
        # create a set of points that span 10 cm (more than enough to span the
        # range of the shelf).
        #mask = gfile["RBBBS"] > gfile["RMAXIS"]
        #f_sep = interp1d(gfile["ZBBBS"][mask], gfile["RBBBS"][mask])
        #romp_sep = f_sep(gfile["ZMAXIS"])
        #romp = np.linspace(romp_sep, romp_sep+0.1, 100)

        # Now see what range of psin these cover.
        #f_psin = RectBivariateSpline(gfile["R"], gfile["Z"], gfile["PSIRZ_NORM"])
        print("Creating psin(R, Z) interpolation function...")
        f_psin = Rbf(Rs, Zs, gfile["PSIRZ_NORM"], epsilon=0.00001)
        #psin_omp = f_psin(romp, np.full(len(romp), gfile["ZMAXIS"]), grid=False)

        # For the shelf see what range of psin it covers.
        shelf_rs = np.linspace(shelf_r1, shelf_r2, 100)
        psin_shelf = f_psin(shelf_rs, np.full(len(shelf_rs), shelf_z))

        # Only need values where psin >= 1.
        mask = psin_shelf >= 1.0
        psin_shelf = psin_shelf[mask]
        shelf_rs = shelf_rs[mask]

        # Calculate the psin range at the X-point.
        romp_range = np.linspace(gfile["RMAXIS"], 2.40, 500)
        rxp_range = np.linspace(gfile["Rx1"], 2.40, 500)
        zomp_range = np.full(len(romp_range), gfile["ZMAXIS"])
        zxp_range = np.full(len(rxp_range), gfile["Zx1"])
        psin_omp = f_psin(romp_range, zomp_range)
        psin_xp = f_psin(rxp_range, zxp_range)
        f_psin_omp = interp1d(psin_omp, romp_range)
        f_psin_xp = interp1d(psin_xp, rxp_range)

        # See what spatial range (R coordinates) that the flux surfaces span
        # at the X-point. Restrict to OMP side.
        #mask = gfile["R"] > gfile["RMAXIS"]
        #print("Generating R(psin, Z) interpolation function...")
        #f_romp = Rbf(gfile["PSIRZ_NORM"][mask], Zs[mask], Rs[mask], epsilon=0.00001)
        romp = f_psin_omp(psin_shelf)
        rxp = f_psin_xp(psin_shelf)

        # The ratio between the spatial distance between the flux surfaces
        # is the flux expansion along the target.
        widths_shelf = [shelf_rs[i+1]-shelf_rs[i] for i in range(0, len(shelf_rs)-1)]
        widths_omp = [romp[i+1]-romp[i] for i in range(0, len(romp)-1)]
        widths_xp = [rxp[i+1]-rxp[i] for i in range(0, len(rxp)-1)]
        #flux_expansion = np.array(widths_shelf) / np.array(widths_omp)
        xp_to_shelf_fe = np.array(widths_xp) / np.array(widths_shelf)
        omp_to_xp_fe = np.array(widths_omp) / np.array(widths_xp)
        rmrs_shelf = shelf_rs - rvsout
        rmrs_xp = rxp_range - gfile["Rx1"]

        fig, ax = plt.subplots()
        ax.plot(rmrs_shelf[:-1], omp_to_xp_fe)
        ax.set_xlabel("R-Rsep (m)", fontsize=14)
        ax.set_ylabel("Flux Expansion, OMP --> XP", fontsize=14)
        fig.tight_layout()
        fig.show()

        return {"r_targ":shelf_rs, "rmrs_targ":rmrs_shelf[:-1], "rmrs_xp":rmrs_xp,
            "xp_to_shelf_fe":xp_to_shelf_fe, "omp_to_xp_fe":omp_to_xp_fe}

    def load_shelf_lp(self, shot, tstart, tend, bins=100):
        """
        Load the LP data just for the data along the shelf.

        Inputs
        shot (int): Shot number.
        tstart (float): Start time in ms for the data.
        tend (float): End time in ms for the data.
        bins (int): Number of bins to bins the data into for cleaning up the data.
        """

        import get_lp
        from scipy.signal import medfilt

        lpdict = get_lp.plot_lps(shot, tstart, tend, bins=bins, tunnel=False, showplot=False)
        labels = ["S-1", "S-2", "S-3", "S-4", "S-5", "S-6", "S-7", "S-8", "S-9",
            "S10", "S11"]

        colors = {}
        for i in range(0, len(labels)):
            colors[labels[i]] = "C{}".format(i)

        data = {"rmrs":[], "ne":[], "te":[], "color":[]}
        for i in range(0, len(lpdict["labels"])):

            label = lpdict["labels"][i].strip()
            if label in labels:
                data["rmrs"].append(lpdict["rminrsep"][i])
                data["ne"].append(lpdict["ne (cm-3)"][i] * 1e6)
                data["te"].append(lpdict["Te (eV)"][i])
                data["color"].append(colors[label])

        sort_idx = np.argsort(data["rmrs"])
        data["rmrs"] = np.array(data["rmrs"])[sort_idx]
        data["ne"] = np.array(data["ne"])[sort_idx]
        data["te"] = np.array(data["te"])[sort_idx]
        data["color"] = np.array(data["color"])[sort_idx]

        # Median filter to the data.
        nefilt = medfilt(data["ne"], 35)
        tefilt = medfilt(data["te"], 35)
        data["nefilt"] = nefilt
        data["tefilt"] = tefilt

        # Peak ne occurs at...
        peak_idx = np.argmax(nefilt)
        rpeak = data["rmrs"][peak_idx]
        nepeak = nefilt[peak_idx]
        print("Peak value at {:.2f} cm from the strike point".format(rpeak*100))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

        ax1.scatter(data["rmrs"], data["ne"], c=data["color"], zorder=5)
        ax1.plot(data["rmrs"], data["nefilt"], color="k", zorder=10)
        ax1.scatter(rpeak, nepeak, marker="*", s=200, color="r", edgecolors="k", zorder=15)
        ax1.set_xlabel("R-Rsep (m)", fontsize=14)
        ax1.set_ylabel("ne (m-3)", fontsize=14)

        ax2.scatter(data["rmrs"], data["te"], c=data["color"])
        ax2.plot(data["rmrs"], data["tefilt"], color="k")
        ax2.set_xlabel("R-Rsep (m)", fontsize=14)
        ax2.set_ylabel("Te (eV)", fontsize=14)

        fig.tight_layout()
        fig.show()

        self.lpdata = data
        return data

    def load_mafot(self, mafot_path):
        """
        Load in a MAFOT run where the connection lengths have been calculated
        for a line across the X-point. These are used to determine how far the
        blob needs to travel before depositing on the target.
        """

        # Load in the mafot data into a dataframe.
        columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
          "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
        df = pd.read_csv(mafot_path, skiprows=52, names=columns, delimiter="\t")
        self.mafot = df
        return df

    def load_xp_r_loc(self, shot, time):
        """
        If mafot is loaded, load the R coordinate of the X-point (it is
        assumed MAFOT was ran along a horizontal line passing through the
        X-point).

        Inputs
        shot (int): Shot number.
        time (float): Time to get location at.
        """

        from gadata import gadata
        import MDSplus
        conn = MDSplus.Connection("atlas.gat.com")
        gaobj = gadata("RXPT1", shot, connection=conn)
        times = gaobj.xdata
        rxpts = gaobj.zdata

        idx = np.argmin(np.abs(times-time))
        rxpt = rxpts[idx]
        print("X-point at R = {:.3} m".format(rxpt))
        self.rxpt = rxpt
        return rxpt

    def run_blp(self, tsdata, lpdata, flux_exp, vr_mean, vr_std, vr_skew,
        xp_loc, nparts, prof_coef=1.0, dist_type="skewnorm", mu_pois=1.0,
        temin=1.0, qtim=1e-7):
        """
        Run the BlobbyLP (blp) Monte Carlo model to intepretively model the
        LP target patterns. A number of assumptions are built into this model.

        - Wherever blobs are generated from, it is assumed that they extend
            down to the X-point.
        - The distance the blob needs to go till hitting the target is thus
            the distance of the X-point.
        - The parallel speed of the blobs is the sound speed.
        - Once a radial velocity is chosen the blobs stay at that velocity.

        Inputs
        tsdata (dict): The dictionary returned from calc_te_params.
        lpdata (dict): The Langmuir probe data as returned from load_shelf_lp.
        flux_exp (dict): The dictionary returned from calc_flux_expansion.
        vr_mean (float): The mean value of the radial blob velocity.
        vr_std (float): The standard deviation of the radial blob velocity distribution.
        vr_skew (float): The skewness of the radial blob velocity distribution.
        xp_loc (float): Parallel distance from the target to the X-point in meters.
        nparts (int): Number of particles to run for the simulation.
        prof_coef (float): Coefficient in front of the resulting target profile.
        dist_type (str): One of "skewnorm" or "poisson". The distribution that
            vr is pulled from.
        mu_pois (float): The mu parameter for the poisson distribution.
        temin (float): Minimum value of Te (bottoms out at this value).
        """

        from tqdm import tqdm
        from scipy.stats import skewnorm, poisson
        from scipy.interpolate import interp1d

        try:
            self.mafot
            mafot_loaded = True
        except:
            print("Warning: Assuming constant connection length of {} m. Please run MAFOT.".format(xp_loc))
            mafot_loaded = False

        if mafot_loaded:
            try:
                self.rxpt
            except:
                print("Error! If using MAFOT, must load in X-point R " + \
                    "coordinate via load_xp_r_loc. Please run that first.")
                sys.exit()

        # Establish inputs to the model.
        tesep = tsdata["tesep"]
        lambda_te = tsdata["lambda_te"] / 100  # cm to m
        rs = np.linspace(0, 0.5, 500)

        # Calculate sound speed. We use the flux expansion between the OMP (where
        # the Te measurements have been mapped to) and the X-point (where we
        # want the measurements at) to expand the lambda_te accordingly. It
        # should be appropriate to just use the average flux expansion since Te
        # measurements have enough uncertainty in them as it is.
        mi = 931.49e6
        avg_omp_to_xp_fe = flux_exp["omp_to_xp_fe"].mean()
        print("Average OMP to X-point flux expansion: {:.2f}".format(1/avg_omp_to_xp_fe))
        te = tesep * np.exp(-rs/(lambda_te/avg_omp_to_xp_fe))
        te = np.clip(te, temin, None)
        cs = np.sqrt((2*te)/mi) * 3e8
        f_cs = interp1d(rs, cs)

        # Using MAFOT, find the length of each field line along each of the
        # R coordinates we will consider.
        rs_at_xp = self.rxpt + rs
        xp_conns = []
        for r in rs_at_xp:
            idx = np.argmin(np.abs(self.mafot["R (m)"] - r))
            xp_conns.append(self.mafot["Lconn (km)"].iloc[idx] * 1000)
        f_xp_conns = interp1d(rs, xp_conns)

        # Choose blobs vr's. Anything below zero remove.
        if dist_type == "skewnorm":
            launch_vrs = skewnorm.rvs(a=vr_skew, loc=vr_mean, scale=vr_std,
                size=nparts)
        elif dist_type == "poisson":
            launch_vrs = poisson.rvs(mu=vr_mean, loc=0, size=nparts)
        below_zero = np.where(launch_vrs <= 0)
        launch_vrs = np.delete(launch_vrs, below_zero)
        print("Actual particles launched: {}".format(len(launch_vrs)))

        targ_locs = np.zeros(len(launch_vrs))
        warn_count = 0
        warn_flag = False
        for i in tqdm(range(0, len(launch_vrs))):

            # Assume all blobs extend to X-point, so that's the distance they need to
            # travel before hitting the target. Follow particle in 2D sense.
            xpos = 0.0
            rpos = 0.0
            vr = launch_vrs[i]
            while True:

                # Perform radial step, see what nearest cs is (which we treat as
                # v parallel).
                rpos = rpos + vr * qtim
                ridx = np.argmin(np.abs(rs-rpos))
                if ridx == len(rs):
                    warn_flag = True
                    warn_count += 1
                vpar = cs[ridx]
                #vpar = f_cs(rpos)

                # Perform parallel step, negative because towards target.
                xpos = xpos + vpar * qtim

                # If xpos is greater than the field line length, we've hit the
                # target. Record location.
                # Since we only calculate the distance from the X-point to the
                # floor, the flux expansion is not that serious and is
                # approximated as 1, i.e., neglected.
                if mafot_loaded:
                    break_condition = xp_conns[ridx]
                    #break_condition = f_xp_conns(rpos)
                else:
                    break_condition = xp_loc
                if xpos > break_condition:
                    targ_locs[i] = rpos
                    break

        if warn_flag:
            print("Warning! Hit end of R range {} times.".format(warn_count))

        # Bin data into histogram, normalize.
        bin_locs = list(np.histogram(targ_locs, 50, density=False))
        bin_locs[0] = bin_locs[0] / bin_locs[0].max()
        plotx = np.array([(bin_locs[1][i]+bin_locs[1][i+1])/2 for i in range(0, len(bin_locs[1])-1)])

        # Also include the vr pdf.
        if dist_type == "skewnorm":
            plot_vr = np.linspace(0, vr_mean*100, 10000)
            pdf = skewnorm.pdf(plot_vr, a=vr_skew, loc=vr_mean, scale=vr_std)
            cutoff = np.where(pdf >= 0.1 * pdf.max())
            plot_vr = plot_vr[cutoff]
            pdf = pdf[cutoff]
        elif dist_type == "poisson":
            #pdf = poisson.pmf(plot_vr, mu=vr_mean, loc=0)
            vr_bins = np.histogram(launch_vrs, 100, density=True)
            plot_vr = np.array([(vr_bins[1][i]+vr_bins[1][i+1])/2 for i in range(0, len(vr_bins[1])-1)])
            pdf = vr_bins[0]


        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

        max_ne = lpdata["nefilt"].max()
        ax1.scatter(lpdata["rmrs"]*100, lpdata["ne"]/max_ne, c=lpdata["color"], zorder=5)
        ax1.plot(lpdata["rmrs"]*100, lpdata["nefilt"]/max_ne, color="k", zorder=10)
        ax1.set_xlabel("R-Rsep (cm)", fontsize=14)
        ax1.set_ylabel("ne (m-3)", fontsize=14)
        #ax1.set_ylabel("Counts", fontsize=14)
        ax1.plot(plotx*100, prof_coef*bin_locs[0], color="r", lw=2, zorder=20)

        ax2.plot(plot_vr, pdf)
        ax2.set_xlabel("vr (m/s)", fontsize=14)
        ax2.set_ylabel("PDF", fontsize=14)

        fig.tight_layout()
        fig.show()
