import sys
import pickle
import openadas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp1d, griddata


class BlobbyFarSOL():

    def __init__(self, shot, plunge):
        self.shot = shot
        self.plunge = plunge

    def load_rcp_profiles(self, path, start_idx=None, end_idx=None, shift=0.0, temod=1.0, timult=1.0):
        print("Loading RCP data...")
        columns = ["Time (ms)", "R (cm)", "Rho", "Isat (A)", "ne (1e18 m-3)", "Te (eV)", "Qpar (MW/m2)", "Mach",
                   "vflow (m/s)", "Vf1 (V)", "Vf2 (V)", "Vf3 (V)"]
        rcp = pd.read_csv(path, names=columns, delimiter="\t", header=0)

        # Some plunges require some additional post-processing, as per Dmitry's instructions. The instructions only
        # apply to the SAS-VW shots. They are:
        # 1. Ignore all points further out than 233 cm, they are dominated by noise in both Te and Isat.
        # 2. If Isat and ne go negative do the following: of 3-4 closest points inward of 233 cm pick the most negative
        #    density and add this offset to the whole profile, so that all density points are positive. The resulting
        #    density will still be on the lower side of the real one (since in the point we set to zero the density is
        #    actually positive), but probably within the error bars.
        if self.shot == 190411:
            if self.plunge == 1:
                rcp["ne (1e18 m-3)"] += 0.724742
            elif self.plunge == 2:
                rcp["ne (1e18 m-3)"] += 0.116934
        elif self.shot == 190440:
            if self.plunge == 1:
                # rcp["ne (1e18 m-3)"] += -rcp["ne (1e18 m-3)"].iloc[8]
                rcp["ne (1e18 m-3)"] += 0.198156
        elif self.shot == 190484:
            if self.plunge == 1:
                # rcp["ne (1e18 m-3)"] += -rcp["ne (1e18 m-3)"].iloc[9]
                rcp["ne (1e18 m-3)"] += 0.719046
            elif self.plunge == 2:
                # rcp["ne (1e18 m-3)"] += -rcp["ne (1e18 m-3)"].iloc[8]
                rcp["ne (1e18 m-3)"] += 0.529617

        if type(end_idx) == type(None):
            rcp = rcp.iloc[start_idx:]
        else:
            rcp = rcp.iloc[start_idx:-end_idx]
        self.rcp = rcp
        self.rcp_r = rcp["R (cm)"].values + shift
        self.rcp_rho = rcp["Rho"].values
        self.rcp_m = rcp["Mach"].values
        self.rcp_te = rcp["Te (eV)"].values * temod
        self.rcp_cs = np.sqrt((self.rcp_te + timult * self.rcp_te) / (2 * 931.49e6)) * 3e8
        self.rcp_ne = rcp["ne (1e18 m-3)"].values * 1e18
        self.rcp_isat = rcp["Isat (A)"].values

    def load_rcp_blobs(self, path, shift=0.0, min_npeaks=0, temod=1.0, timult=1.0, vr_time_avg_mult=1.0):
        print("Loading blob data...")
        columns = ["Time (ms)", "R (cm)", "R-Rsep (cm)", "Isat (A)", "ne (1e18 m-3)", "Te (eV)", "Epol (V/m)",
                   "Erad (V/m)", "vr (m/s)", "Npeaks", "Tblob (1e-6 s)", "Drad (cm)", "Dtot (cm)"]
        blob = pd.read_csv(path, names=columns, delimiter="\t", header=0)
        timestep = (blob["Time (ms)"].iloc[1] - blob["Time (ms)"].iloc[0]) / 1000  # ms to s

        # As above, so below.
        if self.shot == 190440:
            if self.plunge == 1:
                blob["ne (1e18 m-3)"] += 0.198156
        elif self.shot == 190484:
            if self.plunge == 1:
                blob["ne (1e18 m-3)"] += 0.719046
            elif self.plunge == 2:
                blob["ne (1e18 m-3)"] += 0.529617

        # Only accept entries with a minimum (statistically significant) amount of blob counts (e.g., Npeaks).
        blob = blob[blob["Npeaks"] >= min_npeaks]

        # Negative densities indicate bad data, don't include.
        blob = blob[blob["ne (1e18 m-3)"] > 0]

        # This is the time-averaged blob velocity, where we have effectively assumed the blobs to be square waves.
        # self.blob_vr_time_avg = blob["vr (m/s)"].values * blob["Tblob (1e-6 s)"].values * 1e-6 * \
        #     blob["Npeaks"].values / timestep * vr_time_avg_mult

        # Time-averaged blob velocity. We create a synthetic signal of repeated Gaussians, and then take the average
        # (via the integral divided by the period). This process seems to report somewhat lower values for the time
        # averaged vr compared to approximating the signal as a square wave. Ultimately seems to have a negligible
        # effect on the end results though.
        tblob = blob["Tblob (1e-6 s)"].values * 1e-6
        synth_t = np.linspace(0, timestep, 5000)
        time_avg_vr = np.zeros(len(tblob))
        for i in range(0, len(blob["Npeaks"].values)):
            n = int(blob["Npeaks"].values[i])
            vr = blob["vr (m/s)"].values[i]

            # Convert from FWHM to sigma (std. dev.). See same question here:
            # https://stackoverflow.com/questions/10623448/making-gaussians-constrained-by-the-fwhm
            fwhm = tblob[i]
            sigma = fwhm * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

            synth_vr = np.zeros(synth_t.shape)
            for j in range(0, n):
                center_t = (j + 1) * timestep / (n + 1)
                synth_vr += vr * np.exp(-np.square((synth_t - center_t) / sigma))

            # Now compute the mean of the signal.
            time_avg_vr[i] = 1 / timestep * np.trapz(synth_vr, synth_t)
        self.blob_vr_time_avg = time_avg_vr

        self.blob_vr = blob["vr (m/s)"].values
        self.blob_r = blob["R (cm)"].values + shift
        self.blob_te = blob["Te (eV)"].values * temod
        self.blob_ne = blob["ne (1e18 m-3)"].values * 1e18
        self.blob_rmrs = blob["R-Rsep (cm)"].values
        self.blob_drad = blob["Drad (cm)"].values
        self.blob_f = blob["Npeaks"].values / timestep
        self.blob_epol = blob["Epol (V/m)"].values
        self.blob_twidth = blob["Tblob (1e-6 s)"] * 1e-6
        self.blob_npeaks = blob["Npeaks"].values
        self.blob_cs = np.sqrt((self.blob_te + timult * self.blob_te) / (2 * 931.49e6)) * 3e8
        self.blob = blob

    def load_mafot(self, itf_path, otf_path):
        print("Loading MAFOT data...")
        columns = ["R (m)", "Z (m)", "N_toroidal", "Lconn (km)", "psimin",
                   "psimax", "psiav", "pitch angle", "yaw angle", "theta", "psi"]
        mafot_itf = pd.read_csv(itf_path, skiprows=52, names=columns, delimiter="\t")
        mafot_otf = pd.read_csv(otf_path, skiprows=52, names=columns, delimiter="\t")
        self.conns_r = mafot_itf["R (m)"] * 100
        self.conns_l_itf = mafot_itf["Lconn (km)"].values * 1000  # km to m
        self.conns_l_otf = mafot_otf["Lconn (km)"].values * 1000

    def load_gfile(self, path):
        """
        Load the pickled gfile output from the EFIT OMFIT module. We use this to get the magnetic field at the RCP
        location, so we can calculate a more accurate vr (vr = Epol / BT).
        """
        print("Loading gfile...")
        with open(path, "rb") as f:
            self.gfile = pickle.load(f)

    def fix_blob_vr(self):
        """
        The vr's given in the RCP data are assuming a constant BT, using the value from Reviewplus. This leads to an
        underestimation of vr, since BT is smaller out in the far-SOL, so we use the gfile to get a slightly more
        accurate vr (vr = Epol/BT).
        """

        # Interpolate onto the gfile grid to get BT at each blob data point.
        R, Z = np.meshgrid(self.gfile["R"], self.gfile["Z"])
        BT = self.gfile["Bt"].flatten().T
        Bp = self.gfile["Bp"].flatten().T
        gfile_points = np.array((R.flatten(), Z.flatten())).T
        blob_points = np.array((self.blob_r / 100, np.full(len(self.blob_r), -0.185))).T
        self.blob_bt = griddata(gfile_points, BT, blob_points)
        self.blob_bp = griddata(gfile_points, Bp, blob_points)
        self.blob_b = np.sqrt(np.square(self.blob_bt) + np.square(self.blob_bp))

        # Multiply by the factor change. Generally fact > 1 and so vr increases.
        prev_bt = self.blob_epol[0] / self.blob_vr[0]
        fact = prev_bt / self.blob_bt
        self.blob_vr = self.blob_vr * fact
        self.blob_vr_time_avg = self.blob_vr_time_avg * fact
        self.blob_drad = self.blob_vr * self.blob["Tblob (1e-6 s)"] * 1e-6 * 100  # m to cm

    def replace_vr_with_linear(self):
        """
        This swaps out the saved blob_vr_time_avg values with the linear fit values so that experimental noise can
        be dampened some.
        """

        z = np.polyfit(self.blob_r, self.blob_vr_time_avg, 1)
        p = np.poly1d(z)
        self.blob_vr_time_avg = p(self.blob_r)
        if z[0] > 0:
            print("ERROR! blob_vr_time_avg slope is positive!")

    def solve_simple_m(self, window_size=3.0, showplot=True, mach_mult=1.0, conn_run_avg=False,
                       conn_run_avg_window=1.0, neut_mult=1.0):
        """
        Solve Eq. 1.46 in Peter's book for the Mach number. Values from the RCP and MAFOT runs are used as needed.

        mach_mult: Needed to make sure a positive Mach number means towards the OT.
        conn_run_avg (bool): Whether or not to apply a running average to the connection length data.
        """

        # Interpolate the connection lengths onto the RCP domain. Sharp changes in connection length can cause a
        # relatively large change in the scatter of the data, so we provide and option to do a running average
        # of the connection length, essentially smearing it out some and lessing the effect an abrupt change in length
        # can have on the data.
        if conn_run_avg:
            from scipy.ndimage import uniform_filter1d

            # Data is evenly spaced, so let's see how many points covers the specified window.
            step_size = self.conns_r[1] - self.conns_r[0]
            N = int(conn_run_avg_window / step_size)
            f_itf = interp1d(self.conns_r, uniform_filter1d(self.conns_l_itf, N))
            f_otf = interp1d(self.conns_r, uniform_filter1d(self.conns_l_otf, N))
        else:
            f_itf = interp1d(self.conns_r, self.conns_l_itf)
            f_otf = interp1d(self.conns_r, self.conns_l_otf)
        rcp_lconn_itf = f_itf(self.rcp_r)
        rcp_lconn_otf = f_otf(self.rcp_r)

        # Treat the inner target as s=0. This means a positive Mach flow is towards the OTF.
        rcp_s = rcp_lconn_itf
        rcp_smax = rcp_lconn_otf + rcp_lconn_itf

        rcp_lnne = np.log(self.rcp_ne)
        if showplot:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
            ax1.scatter(self.rcp_r, rcp_lnne)

        # Do a running window of lambda_ne at each location.
        lambda_ne = np.zeros(len(self.rcp_r))
        dndr = np.zeros(len(self.rcp_r))
        for i in range(0, len(self.rcp_r)):
            window = np.full(len(self.rcp_r), False)
            r = self.rcp_r[i]
            mask = np.abs(self.rcp_r - r) < window_size / 2
            window[mask] = True
            z = np.polyfit(self.rcp_r[window], rcp_lnne[window], 1)
            p = np.poly1d(z)
            lambda_ne[i] = -1 / z[0]
            if showplot:
                ax1.plot(self.rcp_r[window], p(self.rcp_r[window]), color="k", lw=1)

            # Calculate dndr up front since we have lambda_ne and the fit right here.
            dndr[i] = z[0] * np.exp(z[0] * self.rcp_r[i] + z[1])

        self.dndr = dndr

        # First calculate what the radial blob flux is.
        blob_rad_flux = self.blob_ne * self.blob_vr_time_avg
        self.blob_rad_flux = blob_rad_flux

        # Fit an exponential to extract the decay length and to calculate the radial derivative.
        z_dg = np.polyfit(self.blob_r / 100, np.log(blob_rad_flux), 1)
        dgamma_blob_dr = z_dg[0] * np.exp(z_dg[0] * self.blob_r / 100 + z_dg[1])
        self.dgamma_blob_dr = dgamma_blob_dr

        # Create an interpolation function for a common domain later.
        f_dgamma_blob_dr = interp1d(self.blob_r, dgamma_blob_dr)
        def exp_dgamma_blob_vr(r):

            # We don't really want this extrapolating beyond the measured value range, so cap/floor the values at the
            # last known values.
            # if r < min(self.blob_r):
            #     return z_dg[0] * np.exp(z_dg[0] * min(self.blob_r) / 100 + z_dg[1])
            # elif r > max(self.blob_r):
            #     return z_dg[0] * np.exp(z_dg[0] * max(self.blob_r) / 100 + z_dg[1])
            # else:
            return z_dg[0] * np.exp(z_dg[0] * r / 100 + z_dg[1])

        # Alternatively try a linear fit, less erratic.
        z_dg2 = np.polyfit(self.blob_r / 100, blob_rad_flux, 1)
        def lin_dgamma_blob_vr(r):
            return z_dg2[0]

        # fig, ax = plt.subplots()
        # ax.scatter(self.blob_r / 100, blob_rad_flux)
        # ax.plot(self.blob_r/100, np.exp(z_dg[0] * self.blob_r / 100 + z_dg[1]))
        # ax.plot(self.blob_r/100, z_dg2[0] * self.blob_r / 100 + z_dg2[1])
        # fig.tight_layout()
        # fig.show()

        # A plot of blob vr's to extract a linear fit the blob slope.
        z = np.polyfit(self.blob_r, self.blob_vr_time_avg, 1)
        p = np.poly1d(z)
        dvdr = z[0] * 100
        if showplot:
            ax2.scatter(self.blob_r, self.blob_vr_time_avg)
            ax2.plot(self.blob_r, p(self.blob_r))

        if showplot:
            ax1.set_xlabel("R (cm)")
            ax1.set_ylabel("ln(ne)")
            ax2.set_ylabel("Blob Time-Averaged vr (m/s)")
            fig.tight_layout()
            fig.show()

        def eq(M, s, C, cs):
            return 2 * np.arctan(M) - M - C * s / cs

        # For getting the ionization rate coefficients.
        oa = openadas.OpenADAS()
        rate_df = oa.read_rate_coef_unres("/Users/zamperini/My Drive/Research/Data/openadas/scd96_h.dat")
        rate_df_cx = oa.read_rate_coef_unres("/Users/zamperini/My Drive/Research/Data/openadas/ccd96_d.dat")

        C = (np.pi / 2 - 1) * self.rcp_cs / (rcp_smax / 2)
        # C = np.zeros(len(self.rcp_r))
        vr_solved_noexp_m = C * lambda_ne / 100
        D_solved = np.zeros(len(self.rcp_r))
        C_solved = np.zeros(len(self.rcp_r))
        vr_solved = np.zeros(len(self.rcp_r))
        mach_solved = np.zeros(len(self.rcp_r))
        mach_linear = np.zeros(len(self.rcp_r))
        C_solved2 = np.zeros(len(self.rcp_r))
        vr_solved2 = np.zeros(len(self.rcp_r))
        E_solved = np.zeros(len(self.rcp_r))
        E_solved2 = np.zeros(len(self.rcp_r))
        vr_solved3 = np.zeros(len(self.rcp_r))
        neut_dens5 = np.zeros(len(self.rcp_r))
        neut_dens6 = np.zeros(len(self.rcp_r))
        neut_dens7 = np.zeros(len(self.rcp_r))
        f_cs = interp1d(self.rcp_r, self.rcp_cs)
        f_s = interp1d(self.rcp_r, rcp_s)
        f_smax = interp1d(self.rcp_r, rcp_smax)
        f_lambda_ne = interp1d(self.rcp_r, lambda_ne)
        f_M = interp1d(self.rcp_r, self.rcp_m * mach_mult)
        for i in range(0, len(self.rcp_r)):

            # If greater than zero, we're past the stagnation point starting from IT.
            meas_mach = self.rcp_m[i] * mach_mult
            if meas_mach > 0:
                slope = (meas_mach - 1) / (rcp_s[i] - rcp_smax[i])
            else:
                slope = (meas_mach + 1) / (rcp_s[i] - 0)
            stag_s = -meas_mach / slope + rcp_s[i]
            it_conn = stag_s
            ot_conn = rcp_smax[i] - stag_s
            s = np.abs(rcp_s[i] - stag_s)

            # Other way of calculating s is to just have it how far from the halfway point are we.
            # if rcp_s[i] > rcp_smax[i] / 2:
            #    s2 = rcp_s[i] - rcp_smax[i] / 2
            # else:
            #    s2 = rcp_smax[i] / 2 - rcp_smax[i]

            # Depending on what region we're in, assign the correct value to L and calculate C.
            # if meas_mach > 0:
            #     C[i] = (np.pi / 2 - 1) * self.rcp_cs[i] / ot_conn
            # else:
            #     C[i] = (np.pi / 2 - 1) * self.rcp_cs[i] / it_conn

            # rcp_s = 0 is the inner target.
            M = rcp_s[i] / (rcp_smax[i] / 2) - 1
            mach_linear[i] = rcp_s[i] / (rcp_smax[i] / 2) - 1
            dMds = 1 / (rcp_smax[i] / 2)  # L = smax / 2
            cs = self.rcp_cs[i]
            ioniz_rate_coef = oa.get_rate_coef(rate_df, self.rcp_te[i], self.rcp_ne[i])
            neut_dens5[i] = cs / ioniz_rate_coef * dMds * (1 - (2 * M**2) / (1 + M**2))
            neut_dens6[i] = cs / ioniz_rate_coef * slope * (1 - (2 * meas_mach ** 2) / (1 + meas_mach ** 2))

            # Using the same equation, let's extract a Dperp by using the measured M and backing out what C (and thus D) is.
            # C_solved[i] = (2 * np.arctan(self.rcp_m[i]*mach_mult) - self.rcp_m[i]*mach_mult) * self.rcp_cs[i] / s
            C_solved[i] = (2 * np.arctan(np.abs(self.rcp_m[i])) - np.abs(self.rcp_m[i])) * self.rcp_cs[i] / s
            D_solved[i] = C_solved[i] * np.square(lambda_ne[i] / 100)

            # Similarly solve for a vr, here we assumed dv/dr ~ 0.
            vr_solved[i] = (C_solved[i] + dvdr) * lambda_ne[i] / 100

            # As another approach, we revert back to the simple approximation that the Mach number linearly approaches
            # each target. We allow that the stagnation point not be halfway between the targets. To find the stagnation
            # point, we use the measured M and slope to see at what s value does M = 0. Once we know the stagnation
            # point, we can then do rcp_s - stag_s to estimate how far from the stagnation point the measurement is. That
            # then lets us calculate C from Eq. 1.42 (C = (2arctan(M) - M) * cs / s). Then when we know C, we can
            # finally do vr = lambda_n * C.

            # Solved normally. Effectively all we are doing is modifying the assumption that M=0 halfways between the
            # targets, and instead have it zero wherever we linearly extrapolate it to be.
            C_solved2[i] = (2 * np.arctan(np.abs(self.rcp_m[i])) - np.abs(self.rcp_m[i])) * self.rcp_cs[i] / s
            vr_solved2[i] = (C_solved2[i] + dvdr) * lambda_ne[i] / 100
            # vr_solved2[i] = (C_solved2[i]) * lambda_ne[i] / 100

            # One step further then would be to calculate the parallel electric field. E = (B/Bp)kTe/2eL. For now just
            # approximate B/Bp ~ 5
            if meas_mach > 0:
                E_solved[i] = 5 * self.rcp_te[i] / (2 * ot_conn)
            else:
                E_solved[i] = 5 * self.rcp_te[i] / (2 * it_conn)

            # Using this version of s, calculated another predicted Mach number.
            ans = fsolve(eq, 0.5, args=(s, C[i], self.rcp_cs[i]))[0]
            mach_solved[i] = ans * mach_mult
            # vr_solved3[i] = (C[i] + dvdr) * lambda_ne[i] / 100
            vr_solved3[i] = (C[i]) * lambda_ne[i] / 100

            # Horrible programming, just copying pasting the blob stuff from below and replacing it with the RCP equivalent.
            s = rcp_s[i]
            meas_mach = self.rcp_m[i] * mach_mult
            it_conn = stag_s
            ot_conn = rcp_smax[i] - stag_s
            if meas_mach > 0:
                slope = (meas_mach - 1) / (s - rcp_smax[i])
                L = ot_conn
            else:
                slope = (meas_mach + 1) / (s - 0)
                L = it_conn

            tmp_M = self.rcp_m[i] * mach_mult
            # tmp_M = rcp_s[i] / (rcp_smax[i] / 2) - 1  # M is solved above ignoring the measurement, assuming linear increase between targets.
            # tmp_M = slope * s + 1
            # dMds = 1 / (rcp_smax[i] / 2)  # L = smax / 2
            dMds = slope
            n0 = self.rcp_ne[i] * (1 + tmp_M ** 2)
            # dnds = n0 * 2 * L ** 2 * (L - s) / (2 * L ** 2 - 2 * L * s + s ** 2) ** 2
            dnds = -2 * tmp_M / (1 + tmp_M**2)**2 * dMds * n0
            dnds_div_n = -2 * tmp_M / (1 + tmp_M**2) * dMds  # We've plugged in ns0/n = 1 + M^2.
            vpar = tmp_M * self.rcp_cs[i]
            #dvds = self.rcp_cs[i] / L
            dvds = self.rcp_cs[i] * dMds
            print("-----")
            print("r, cs, ne, M, dM/ds, s, smax= {:.2f}  {:.2f}  {:.2e}  {:.2f}  {:.2f}  {:.2f}  {:.2f}".format(self.rcp_r[i], self.rcp_cs[i], self.rcp_ne[i], tmp_M, dMds, rcp_s[i], rcp_smax[i]))
            print("1 / ioniz_coff = {:.2e}".format(1 / ioniz_rate_coef))
            print("dnds * vpar / ne = {}".format(dnds_div_n * vpar))
            print("dvds = {}".format(dvds))
            print("exp_dg / ne = {}".format(exp_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i]))
            print("sum = {:.2e}".format(dnds_div_n * vpar + dvds + exp_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i]))
            #print("lin_dg / ne = {}".format(lin_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i]))

            #neut_dens7[i] = 1 / ioniz_rate_coef * (dnds * vpar / self.rcp_ne[i] + dvds + exp_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i])
            #neut_dens7[i] = 1 / ioniz_rate_coef * (dnds * vpar / self.rcp_ne[i] + dvds + lin_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i])
            neut_dens7[i] = 1 / ioniz_rate_coef * (dnds_div_n * vpar + dvds + exp_dgamma_blob_vr(self.rcp_r[i]) / self.rcp_ne[i])
            #print("{:.2e} x {:.2e} = {:.2e}".format(1 / (ioniz_rate_coef * self.rcp_ne[i]), neut_dens7[i]))

        # Needed to know what range of R is fine.
        min_r = max(self.rcp_r.min(), self.blob_r.min())
        max_r = min(self.rcp_r.max(), self.blob_r.max())
        blob_mask = np.logical_and(self.blob_r > min_r, self.blob_r < max_r)

        # Load LLAMA data if available.
        if self.shot in [190484, 190485, 190486]:

            # 190440 and 190442 not like the others... ignoring for now.
            if self.shot in [190440, 190442]:
                llama_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/LLAMA_190422.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190440/190440_3000.pickle"
            if self.shot == 190484:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"
                a = 8.20e15
                b = 9.59
                c = -5.92e15
            if self.shot == 190485:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190485_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190485/190485_3000.pickle"
                a = -5.46e16
                b = -1.7
                c = 5.67e16
            if self.shot == 190486:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190486_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190486/190486_3000.pickle"
                a = 1.54e16
                b = 4.98
                c = -1.29e16
            llama = np.load(llama_path)
            with open(gfile_path, "rb") as f:
                gfile = pickle.load(f)

            # Get some gfile stuff so we can go from R, Z to psin.
            R = gfile["R"]
            Z = gfile["Z"]
            Rs, Zs = np.meshgrid(R, Z)
            psin = gfile["PSIRZ_NORM"]
            blob_coord = zip(self.blob_r / 100, np.full(len(self.blob_r), -0.185))
            blob_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(blob_coord))

            lrho = np.square(llama["psi_n"])
            lneut = llama["nDens_LFS"]
            lneut_err = llama["nDens_LFS_err"]
            lion = llama["ion_LFS"]
            lion_err = llama["ion_LFS_err"]
            lpsin = llama["psi_n"]

            # Interpolate to determine what the neutral density is on the psin of the RCP measurement.
            # f_neut = interp1d(lpsin, lneut)
            # blob_neut = f_neut(blob_psin)

            # In a separate script we have done exponential fits to the data so that we can extrapolate some. Fit
            # was for values of psin - 1.
            blob_neut = (a * np.exp(b * (blob_psin - 1)) + c) * neut_mult


        # Take the time-averaged blob vr to predict lambda_ne.
        solved_lambda_ne = np.zeros(len(self.blob_r))
        solved_lambda_ne2 = np.zeros(len(self.blob_r))
        solved_lambda_ne3 = np.zeros(len(self.blob_r))
        solved_lambda_ne4 = np.zeros(len(self.blob_r))
        neut_dens2 = np.zeros(len(self.blob_r))
        neut_dens3 = np.zeros(len(self.blob_r))
        neut_dens4 = np.zeros(len(self.blob_r))
        siz = np.zeros(len(self.blob_r))

        ioniz_rate_coeffs = np.zeros(len(self.blob_r))
        cx_rate_coeffs = np.zeros(len(self.blob_r))
        for i in range(0, len(self.blob_r)):

            if blob_mask[i] == False:
                continue
            r = self.blob_r[i]

            # Calculate C, and then using the blob values calculate lambda_ne.
            tmp_C = (np.pi / 2 - 1) * f_cs(r) / (f_smax(r) / 2)
            solved_lambda_ne[i] = self.blob_vr_time_avg[i] / tmp_C * 100
            # solved_lambda_ne[i] = self.blob_vr_time_avg[i] / (tmp_C + dvdr) * 100

            # For this next one, we've derived some equations starting with the 2D conservation equations and pulled
            # a couple of things from the simple SOL.
            # L = f_smax(r) / 2
            s = f_s(r)
            # if f_s(r) > f_smax(r) / 2:
            #     s = f_s(r) - f_smax(r) / 2
            # else:
            #     s = f_smax(r) / 2 - f_s(r)

            # If greater than zero, we're past the stagnation point starting from IT.
            meas_mach = f_M(r)
            it_conn = stag_s
            ot_conn = f_smax(r) - stag_s
            if meas_mach > 0:
                slope = (meas_mach - 1) / (s - f_smax(r))
                L = ot_conn
            else:
                slope = (meas_mach + 1) / (s - 0)
                L = it_conn
            stag_s = -meas_mach / slope + s
            stmp = np.abs(s - stag_s)

            # print("rmrs = {:.2f}: smax/2 = {:.1f}   and  {:.1f}".format(self.blob_rmrs[i], f_smax(r)/2, L))

            # tmp_M = s / L - 1
            tmp_M = f_M(r)
            n0 = self.blob_ne[i] * (1 + tmp_M ** 2)
            dnds = n0 * 2 * L ** 2 * (L - s) / (2 * L ** 2 - 2 * L * s + s ** 2) ** 2
            vpar = tmp_M * self.blob_cs[i]
            dvds = self.blob_cs[i] / L

            # First calculate what lambda would be without any ionization.
            tmp = dnds * vpar / (self.blob_ne[i] * self.blob_vr_time_avg[i]) + 1 / self.blob_vr_time_avg[i] * dvds \
                  + 1 / self.blob_vr_time_avg[i] * dvdr
            solved_lambda_ne2[i] = 1 / tmp * 100  # m to cm

            # Now using everything, calculate what the neutral density should be.
            tmp = vpar / self.blob_ne[i] * dnds + dvds + dvdr - self.blob_vr_time_avg[i] / (f_lambda_ne(r) / 100)
            ioniz_rate_coef = oa.get_rate_coef(rate_df, self.blob_te[i], self.blob_ne[i])
            ioniz_rate_coeffs[i] = ioniz_rate_coef
            neut_dens2[i] = tmp / ioniz_rate_coef



            # This derivation had vr and dv/dr cancel out. Suspicious...
            #neut_dens4[i] = self.blob_cs[i] / (L * ioniz_rate_coef) * (dnds * (s - L) / self.blob_ne[i] + 1)
            neut_dens4[i] = self.blob_cs[i] / (L * ioniz_rate_coef) * (1 - (2 * (L-s)**2) / (2*L**2 - 2 * L * s + s**2))

            # These aren't used here, but we calculate them since they're used later in the plotting portion.
            cx_rate_coef = oa.get_rate_coef(rate_df_cx, self.blob_te[i], self.blob_ne[i])
            cx_rate_coeffs[i] = cx_rate_coef

            # Alternatively, we can avoid the ADAS rates and just calculate the ionization rate, Siz.
            # print("term1: {:.2e}".format(dnds * vpar))
            # print("term2: {:.2e}".format(dvds * self.blob_ne[i]))
            # print("term3: {:.2e}".format(dvdr * self.blob_ne[i]))
            # print("term4: {:.2e}".format(-self.blob_ne[i] * self.blob_vr_time_avg[i]/ (f_lambda_ne(r) / 100)))
            siz[i] = dnds * vpar + dvds * self.blob_ne[i] + dvdr * self.blob_ne[i] - self.blob_ne[i] \
                     * self.blob_vr_time_avg[i] / (f_lambda_ne(r) / 100)

            # # Newer equations that remove the simple SOL assumption for ne and instead use the momentum equation to
            # # # derive it.
            # # In these, n0 is the RCP value at its corresponding s value. C is a constant that we get by assuming
            # # isothermal flux surfaces.
            # n0 = self.blob_ne[i]
            md = 2 * 931.49e6 / 3e8 ** 2
            # C = md / self.blob_te[i] * dvds - md * self.blob_cs[i] * ioniz_rate_coef / 2
            # print(C)
            #
            # # dn/ds is -C n0 exp(C(s0-s)), where s0 is the s of the RCP measurement. We are evaluating the derivative
            # # at s0, so it is just -C * n0. The equations that follow are unchanged from the above.
            # dnds = -C * n0
            # tmp = dnds * vpar / (self.blob_ne[i] * self.blob_vr_time_avg[i]) + 1 / self.blob_vr_time_avg[i] * dvds \
            #       + 1 / self.blob_vr_time_avg[i] * dvdr
            # solved_lambda_ne3[i] = 1 / tmp * 100
            # tmp = vpar / self.blob_ne[i] * dnds + dvds + dvdr - self.blob_vr_time_avg[i] / (f_lambda_ne(r) / 100)
            # neut_dens3[i] = tmp / ioniz_rate_coef

            # Derived using the momentum equation (not sure this is working correctly).
            tmp = md / self.blob_te[i] * dvdr + 0  # 0 is assuming no neutrals.
            solved_lambda_ne3[i] = 1 / tmp * 100
            tmp = 1 / f_lambda_ne(r) - md / self.blob_te[i] * dvdr
            neut_dens3[i] = tmp / (md * self.blob_vr_time_avg[i] * ioniz_rate_coef)

            # Recalculate as we did for lambda_ne2, except include the neutral term from LLAMA data if it was
            # available.
            if self.shot in [190484, 190485, 190486]:
                tmp = dnds * vpar / (self.blob_ne[i] * self.blob_vr_time_avg[i]) + 1 / self.blob_vr_time_avg[i] * dvds \
                      + 1 / self.blob_vr_time_avg[i] * dvdr - blob_neut[i] * ioniz_rate_coef / self.blob_vr_time_avg[i]
                solved_lambda_ne4[i] = 1 / tmp * 100  # m to cm


        # Using all our experimental data, use it to predict what the neutral density must be.
        f_lambda_ne = interp1d(self.rcp_r, lambda_ne)

        neut_dens = np.zeros(len(self.blob_r))
        for i in range(0, len(self.blob_r)):
            if blob_mask[i] == False:
                continue
            r = self.blob_r[i]
            te = self.blob_te[i]
            ne = self.blob_ne[i]
            ioniz_rate_coef = oa.get_rate_coef(rate_df, te, ne)
            neut_dens[i] = ((np.pi / 2 - 1) * f_cs(r) / (f_smax(r) / 2) - self.blob_vr_time_avg[i] /
                            (f_lambda_ne(r) / 100)) / ioniz_rate_coef

        # Store results.
        self.mach_solved = mach_solved
        self.C_solved = C_solved
        self.D_solved = D_solved
        self.vr_solved = vr_solved
        self.rcp_s = rcp_s
        self.rcp_smax = rcp_smax
        self.lambda_ne = lambda_ne
        self.vr_solved_noexp_m = vr_solved_noexp_m
        self.C_solved2 = C_solved2
        self.vr_solved2 = vr_solved2
        self.dvdr = dvdr
        self.E_solved = E_solved
        self.vr_solved3 = vr_solved3
        self.C = C
        self.solved_lambda_ne = solved_lambda_ne
        self.neut_dens = neut_dens
        self.solved_lambda_ne2 = solved_lambda_ne2
        self.solved_lambda_ne3 = solved_lambda_ne3
        self.neut_dens2 = neut_dens2
        self.siz = siz
        self.blob_ioniz_rate_coeffs = ioniz_rate_coeffs
        self.blob_cx_rate_coeffs = cx_rate_coeffs
        self.solved_lambda_ne4 = solved_lambda_ne4
        self.neut_dens4 = neut_dens4
        self.neut_dens5 = neut_dens5
        self.neut_dens6 = neut_dens6
        self.neut_dens7 = neut_dens7
        self.mach_linear = mach_linear

    def calc_neut(self, verbal=False):
        """
        The previous function had gotten pretty messy that I've decided to start afresh.
        """

        # First calculate what the radial blob flux is.
        blob_rad_flux = self.blob_ne * self.blob_vr_time_avg
        self.blob_rad_flux = blob_rad_flux

        # Fit an exponential to extract the decay length and to calculate the radial derivative.
        z = np.polyfit(self.blob_r / 100, np.log(blob_rad_flux), 1)
        dgamma_blob_dr = z[0] * np.exp(z[0] * self.blob_r / 100 + z[1])
        self.dgamma_blob_dr = dgamma_blob_dr

        # Create an interpolation function for a common domain later.
        f_dgamma_blob_dr = interp1d(self.blob_r, dgamma_blob_dr)

        if verbal:
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.scatter(self.blob_r, blob_rad_flux)
            ax.plot(self.blob_r, np.exp(z[0] * self.blob_r / 100 + z[1]))
            ax.set_xlabel("R (m)")
            ax.set_ylabel("blob_rad_flux (m-2 s-1)")
            fig.tight_layout()
            fig.show()

        # Now the normal RCP data, which is the "left hand side" of our equation.
        # First dndr.
        z = np.polyfit(self.rcp_r / 100, np.log(self.rcp_ne), 1)
        dndr = z[0] * np.exp(z[0] * self.rcp_r + z[1])



    def plot_results(self, summary_range=(None, None), mach_mult=1.0):

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 7), sharex=True)

        ax1.axhline(0, color="k", linestyle="--")
        ax1.plot(self.rcp_r, self.mach_solved, label="Simple SOL")
        ax1.plot(self.rcp_r, self.mach_linear, label="Linear")
        ax1.plot(self.rcp_r, self.rcp_m * mach_mult, label="Measured")
        ax1.set_xlim([224, 235])
        ax1.legend()
        ax1.set_ylabel("Mach")

        ax2.plot(self.conns_r, self.conns_l_itf, label="ITF")
        ax2.plot(self.conns_r, self.conns_l_otf, label="OTF")
        ax2.legend()
        ax2.set_ylabel("Connection Length (m)")
        ax2.set_yscale("log")
        ax2.set_ylim([3, 100])

        ax3.plot(self.rcp_r, self.D_solved, color="k")
        ax33 = ax3.twinx()
        ax33.plot(self.rcp_r, self.vr_solved, color="r", label="Simple SOL")
        ax33.plot(self.blob_r, self.blob_vr_time_avg, color="r", linestyle="--", label="Time-averaged")
        ax33.plot(self.rcp_r, self.vr_solved_noexp_m, color="r", linestyle=":", label="Simple SOL - No Exp M")
        ax3.set_ylabel("Dperp (m2/s)")
        ax33.set_ylabel("vr (m/s)", color="r")
        ax33.tick_params(axis="y", labelcolor="r")
        ax33.legend()

        fig.tight_layout()
        fig.show()

        # An additional plot for the SAS-VW plunges to consider LLAMA data as an order of magnitude comparison.
        # Get LFS neutral density. Aaron said the HFS isn't good data.
        if self.shot in [190440, 190442, 190484, 190485, 190486]:
            if self.shot in [190440, 190442]:
                llama_path = "/Users/zamperini/Documents/d3d_work/divimp_files/190423/LLAMA_190422.npz"
            if self.shot == 190484:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190484_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle"
            if self.shot == 190485:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190485_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190485/190485_3000.pickle"
            if self.shot == 190486:
                llama_path = "/Users/zamperini/Documents/d3d_work/files/LLAMA_190486_.npz"
                gfile_path = "/Users/zamperini/Documents/d3d_work/mafot_files/190486/190486_3000.pickle"
            llama = np.load(llama_path)
            with open(gfile_path, "rb") as f:
                gfile = pickle.load(f)

            # Get some gfile stuff so we can go from R, Z to psin.
            R = gfile["R"]
            Z = gfile["Z"]
            Rs, Zs = np.meshgrid(R, Z)
            psin = gfile["PSIRZ_NORM"]
            rcp_coord = zip(self.blob_r/100, np.full(len(self.blob_r), -0.185))
            rcp_psin = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(rcp_coord))

            # Some files give rho, others psin.
            # if shot in [190485]:
            #     lrho = llama["rho"]
            #     lneut = llama["nDens_LFS"].mean(axis=0)
            #     # lneut_err = np.sqrt(np.square(llama["nDens_LFS_err"]).sum(axis=0))
            #     lneut_err = np.nanmean(llama["nDens_LFS_err"], axis=0)
            #     lion = llama["ion_LFS"].mean(axis=0)
            #     lion_err = np.nanmean(llama["ion_LFS_err"], axis=0)
            # else:
            lrho = np.square(llama["psi_n"])
            lneut = llama["nDens_LFS"]
            lneut_err = llama["nDens_LFS_err"]
            lion = llama["ion_LFS"]
            lion_err = llama["ion_LFS_err"]
            lpsin = llama["psi_n"]

            # Use the normal RCP data to get a function of rho(R).
            min_r = max(self.rcp_r.min(), self.blob_r.min())
            max_r = min(self.rcp_r.max(), self.blob_r.max())
            blob_mask = np.logical_and(np.logical_and(self.blob_r > min_r, self.blob_r < max_r), self.neut_dens2 > 0)
            # f_rho = interp1d(self.rcp_r, self.rcp_rho)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharex=True)
            ax1.axvline(1.0, color="k", zorder=5, linestyle="--")
            # ax1.fill_between(lrho, lneut - lneut_err, lneut + lneut_err, color="tab:pink", alpha=0.3)
            # ax1.errorbar(lpsin, lneut, lneut_err, lw=0, marker=".", ms=10, mfc="tab:pink", mec="k", zorder=15,
            #              elinewidth=1, ecolor="k", label="LLAMA")
            ax1.scatter(lpsin, lneut, marker=".", s=10, color="tab:pink", edgecolors="k", zorder=15, label="LLAMA")
            # ax1.plot(lrho, lneut, label="LLAMA", color="tab:pink", lw=3)
            # ax1.plot(f_rho(self.blob_r[blob_mask]), self.neut_dens2[blob_mask], label="RCP-Derived", color="tab:cyan",
            #          lw=3, zorder=25)
            ax1.scatter(rcp_psin[blob_mask], self.neut_dens2[blob_mask], label="RCP-Derived",
                        color="tab:cyan", zorder=25, s=125, edgecolors="k", marker="*")
            ax1.set_yscale("log")
            ax1.set_ylim([1e14, 1e18])
            ax1.set_xlim([0.99, 1.18])
            ax1.grid(which="both", alpha=0.3)
            ax1.set_ylabel("Neutral Density (m-3)", fontsize=16)
            ax1.set_xlabel("Psin", fontsize=16)
            ax1.tick_params(axis="both", which="major", labelsize=14)
            ax1.legend(fontsize=14, framealpha=1.0).set_zorder(35)

            ax2.axvline(1.0, color="k", zorder=5, linestyle="--")
            # ax2.fill_between(lrho, lion - lion_err, lion + lion_err, color="tab:pink", alpha=0.3)
            # ax2.errorbar(lpsin, lion, lion_err, lw=0, marker=".", ms=10, mfc="tab:pink", mec="k", zorder=15,
            #              elinewidth=1, ecolor="k", label="LLAMA")
            ax2.scatter(lpsin, lion, marker=".", s=10, color="tab:pink", edgecolors="k", zorder=15, label="LLAMA")
            # ax2.plot(lrho, lion, label="LLAMA", color="tab:pink", lw=3)
            # ax2.plot(f_rho(self.blob_r[blob_mask]), self.siz[blob_mask], label="RCP-Derived", color="tab:cyan",
            #          lw=3)
            ax2.scatter(rcp_psin[blob_mask], self.siz[blob_mask], label="RCP-Derived",
                        color="tab:cyan", zorder=25, s=125, edgecolors="k", marker="*")
            ax2.set_yscale("log")
            ax2.set_xlim([0.99, 1.18])
            ax2.set_ylim([1e19, 1e22])
            ax2.grid(which="both", alpha=0.5)
            ax2.set_ylabel("Ionization Rate (m-3 s-1)", fontsize=16)
            ax2.set_xlabel("Psin", fontsize=16)
            ax2.tick_params(axis="both", which="major", labelsize=14)
            ax2.legend(fontsize=14, framealpha=1.0).set_zorder(35)
            fig.tight_layout()
            fig.show()

            # The midplane collisionality.
            # The ei-collision frequency calculation. From page 197 in Friedberg and Eq. 20 in Myra PoP 2006.
            boltz = 8.617e-5  # eV/K
            elec = 1.609e-19
            md = 2.014 * 1.66e-27  # kg
            me = 5.486e-4 * 1.66e-27  # kg
            ln_lambda = 15
            eps = 8.854e-12  # F/m
            f_conn_otf = interp1d(self.conns_r, self.conns_l_otf)
            nu_ei_hat = 1.33e5 * self.blob_ne[blob_mask] * 1e-20 / np.power(self.blob_te[blob_mask] * 1e-3, 3 / 2)
            nu_ei = 1.4 * (5.486e-4 / 2.014) * nu_ei_hat  # 1.4 * (me/mi) * nu_hat
            vi = np.sqrt(3 * elec * self.blob_te[blob_mask] / md)  # m/s
            vte = np.sqrt(2 * elec * self.blob_te[blob_mask] / md)  # m/s
            nu_ei = (1 / (4 * np.pi) * self.blob_ne[blob_mask] * elec ** 4 / (eps ** 2 * me * md) * ln_lambda) * \
                    (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))  # 1/s
            elec_gyrofreq = 1.609e-19 * self.blob_b[blob_mask] / (9.109e-31)
            ion_gyrofreq = 1.609e-19 * self.blob_b[blob_mask] / (2.014 * 1.66e-27)
            ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * self.blob_te[blob_mask]) / (
                    elec * self.blob_b[blob_mask])
            norm_coll = nu_ei * f_conn_otf(self.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.plot(rcp_psin[blob_mask], norm_coll)
            ax.set_xlabel("Psin", fontsize=14)
            ax.set_ylabel("Midplane collisionality", fontsize=14)
            fig.tight_layout()
            fig.show()

        # Print out some info.
        if type(summary_range[0]) == type(None):
            summary_range[0] = self.rcp_r.min()
        if type(summary_range[1]) == type(None):
            summary_range[1] = self.rcp_r.max()
        print("Summary between {}-{} cm".format(summary_range[0], summary_range[1]))
        rcp_mask = np.logical_and(self.rcp_r > summary_range[0], self.rcp_r < summary_range[1])
        blob_mask = np.logical_and(self.blob_r > summary_range[0], self.blob_r < summary_range[1])
        print(" Avg. solved vr:          {:.2f} m/s".format(self.vr_solved[rcp_mask].mean()))
        print(" Avg. solved vr (no exp): {:.2f} m/s".format(self.vr_solved_noexp_m[rcp_mask].mean()))
        print(" Avg. time-averaged vr:   {:.2f}".format(self.blob_vr_time_avg[blob_mask].mean()))
        print(" Avg. Dperp:              {:.2f}".format(self.D_solved[rcp_mask].mean()))


def main(shot, plunge, showplot=True, min_npeaks=5, temod=1.0, usage=1, verbal=False):
    bfs = BlobbyFarSOL(shot, plunge)

    # Some options.
    vr_time_avg_mult = 1.0
    timult = 1.0
    window_size = 4.0
    conn_run_avg = True
    conn_run_avg_window = 0.5

    # Case by case selection of the needed files.
    # USN Rev. BT: +M = Towards OTF
    # USN Fwd. BT: +M = Towards ITF
    # LSN Rev. BT: +M = Towards ITF
    # LSM Fwd. BT: +M = Towards OTF
    if shot == 190411:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190411_1.tab", 0, 5,
                                  temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190411_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_-1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190411/190411_3000.pickle")
            bfs.nesep = 0.0
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/all_plunges/MP190411_2.tab", 15, 1,
                                  temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/all_ca/CA_190411_2.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_-1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190411/lam_rcp3000_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190411/190411_3000.pickle")
            bfs.nesep = 0.0

    if shot == 190440:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190440_1.tab", 12, 3, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190440_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp2000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp2000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190440/190440_2000.pickle")
            bfs.nesep = 5.58e18
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190440_2.tab", 0, 3, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190440_2.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp3000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190440/lam_rcp3000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190440/190440_3000.pickle")
            bfs.nesep = 7.6e18

    elif shot == 190442:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190442_1.tab", 6, 7, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190442_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190442/lam_rcp3000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190442/lam_rcp3000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190442/190442_3000.pickle")
            bfs.nesep = 5.25e18
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190442_2.tab", 6, 7, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190442_2.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190442/lam_rcp4500_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190442/lam_rcp4500_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190442/190442_4500.pickle")
            bfs.nesep = 1.13e19

    elif shot == 190484:
        # shift = 0.7
        shift = 0.0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_1.tab", 20, 6,
                                  temod=temod, shift=shift)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190484_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod, shift=shift)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp3000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp3000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_3000.pickle")
            bfs.nesep = 8.4e18
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190484_2.tab", 20, 2,
                                  temod=temod, shift=shift)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190484_2.tab",
                               min_npeaks=min_npeaks, shift=shift, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp4000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190484/lam_rcp4000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190484/190484_4000.pickle")
            bfs.nesep = 9.28e18

    elif shot == 190485:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_1.tab", 12, 9,
                                  shift=1.0, timult=timult, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190485_1.tab",
                               min_npeaks=min_npeaks, shift=1.0, vr_time_avg_mult=vr_time_avg_mult, timult=timult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190485/lam_rcp3000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190485/lam_rcp3000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190485/190485_3000.pickle")
            bfs.nesep = 1.07e19
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190485_2.tab", 6, 3,
                                  shift=-1.0, timult=timult, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190485_2.tab",
                               min_npeaks=min_npeaks, shift=-1.0, vr_time_avg_mult=vr_time_avg_mult, timult=timult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190485/lam_rcp3980_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190485/lam_rcp3980_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190485/190485_3980.pickle")
            bfs.nesep = 1.07e19

    elif shot == 190486:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190486_1.tab", 10, 2,
                                  shift=1.0, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190486_1.tab",
                               min_npeaks=min_npeaks, shift=1.0, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190486/lam_rcp3000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190486/lam_rcp3000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190486/190486_3000.pickle")
            bfs.nesep = 1.16e19
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/MP190486_2.tab", 8, 1,
                                  shift=-0.5, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/2022-36-03/CA_190486_2.tab",
                               min_npeaks=min_npeaks, shift=-0.5, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/190486/lam_rcp4000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/190486/lam_rcp4000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/190486/190486_4000.pickle")
            bfs.nesep = 1.17e19

    elif shot == 184267:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP184267_1.tab", 5, 2, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_184267_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/184267/lam_rcp1600_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/184267/lam_rcp1600_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/184267/184267_1600.pickle")
            bfs.nesep = 7.47e18
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP184267_2.tab", 10, 2, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_184267_2.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/184267/lam_rcp5000_conn_+1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/184267/lam_rcp5000_conn_-1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/184267/184267_5000.pickle")
            bfs.nesep = 7.39e18

    elif shot == 167193:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP167193_1.tab", 27, 18, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_167193_1.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/167193/lam_rcp2000_conn_-1.dat",
                           # Took care of this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/167193/lam_rcp2000_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/167193/167193_2000.pickle")
            bfs.nesep = 1.15e19
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP167193_2.tab", 26, 9, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_167193_2.tab",
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/167193/lam_rcp3000_conn_-1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/167193/lam_rcp3000_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/167193/167193_3000.pickle")
            bfs.nesep = 1.13e19

    elif shot == 167195:
        shift = -0.015
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP167195_1.tab", 1, 6, shift, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_167195_1.tab", shift,
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp2000_conn_-1.dat",
                           # Took care of this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp2000_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/167195/167195_2000.pickle")
            bfs.nesep = 1.17e19
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP167195_2.tab", 7, 5, shift=-0.5, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_167195_2.tab",
                               shift=-0.5, min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_-1.dat",
                           "/Users/zamperini/Documents/d3d_work/mafot_files/167195/lam_rcp4500_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/167195/167195_4500.pickle")
            bfs.nesep = 1.17e19

    elif shot == 184527:
        shift = 0
        if plunge == 1:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP184527_1.tab", 12, 5, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_184527_1.tab", shift,
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/184527/lam_rcp1600_conn_-1.dat",  # Got this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/184527/lam_rcp1600_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/184527/184527_1600.pickle")
            bfs.nesep = 6.81e18
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP184527_2.tab", 14, 3, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_184527_2.tab", shift,
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/184527/lam_rcp3400_conn_-1.dat",  # Got this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/184527/lam_rcp3400_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/184527/184527_3400.pickle")
            bfs.nesep = 6.42e18

    elif shot == 187111:
        shift = -1.0
        if plunge == 1:  # This one is rough.
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP187111_1.tab", 5, 3, shift=shift, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_187111_1.tab", shift,
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/187111/lam_rcp1600_conn_-1.dat",  # Got this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/187111/lam_rcp1600_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/187111/187111_1600.pickle")
            bfs.nesep = 1.06e19
        elif plunge == 2:
            bfs.load_rcp_profiles("/Users/zamperini/My Drive/Research/Data/rcp_data/MP187111_2.tab", shift=shift, temod=temod)
            bfs.load_rcp_blobs("/Users/zamperini/My Drive/Research/Data/rcp_data/blob_analyzed/CA_187111_2.tab", shift,
                               min_npeaks=min_npeaks, vr_time_avg_mult=vr_time_avg_mult, temod=temod)
            bfs.load_mafot("/Users/zamperini/Documents/d3d_work/mafot_files/187111/lam_rcp3400_conn_-1.dat",  # Got this
                           "/Users/zamperini/Documents/d3d_work/mafot_files/187111/lam_rcp3400_conn_+1.dat")
            bfs.load_gfile("/Users/zamperini/Documents/d3d_work/mafot_files/187111/187111_3400.pickle")
            bfs.nesep = 1.04e19

    if shot in [184527, 187111, 167193, 167195, 190411]:
        mach_mult = -1.0
    else:
        mach_mult = 1.0

    # Apply vr fix first, then do all the solving.
    bfs.fix_blob_vr()
    if usage == 1:
        bfs.replace_vr_with_linear()
        bfs.solve_simple_m(showplot=showplot, window_size=window_size, mach_mult=mach_mult, conn_run_avg=conn_run_avg,
                       conn_run_avg_window=conn_run_avg_window)
    elif usage == 2:
        bfs.calc_neut(verbal=verbal)

    if usage == 1 and showplot:
        if shot in [187111]:
            pass
        else:
            bfs.plot_results(summary_range=[225, 230], mach_mult=mach_mult)
    return bfs


def main_all(smooth_window=None):
    # if type(smooth_window) != type(None):
    from scipy.signal import savgol_filter

    plot_opt = 2
    plot_data = "lambda_ne_neut"  # category, lambda_ne, lambda_ne_neut, neut_dens, vr, myra or instability.

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    fig, ax1 = plt.subplots(figsize=(6, 5))
    x = np.array([])
    y = np.array([])
    r = np.array([])
    te = np.array([])
    ne = np.array([])
    conn = np.array([])
    lambda_ne = np.array([])
    dvdr = np.array([])
    drad = np.array([])
    y_ExB = np.array([])
    M = np.array([])
    vr = np.array([])
    f = np.array([])
    nesep = np.array([])
    twidth = np.array([])
    isat = np.array([])
    opt2x = np.array([])
    opt2y = np.array([])
    opt2c = np.array([])
    opt2m = np.array([])
    ahats = np.array([])
    count = 0
    boltz = 8.617e-5  # eV/K
    elec = 1.602e-19
    md = 2.014 * 1.66e-27  # kg
    me = 5.486e-4 * 1.66e-27  # kg
    ln_lambda = 15
    eps = 8.854e-12  # F/m
    shots = [190411, 190440, 190442, 190484, 190485, 190486, 184267, 167193, 167195, 184527, 187111]
    # shots = [190440, 190442, 190484, 190485, 190486, 184267, 167193, 167195]
    labels = np.array([])
    for shot in shots:
        for plunge in [1, 2]:

            if shot == 187111 and plunge == 1:
                continue
            if shot == 167195 and plunge == 1:
                continue

            # Run normally, then pull out the solved vr and time-averaged blob vr. Grab some other things while we're
            # # here as well.
            print("#{} Plunge #{}".format(shot, plunge))
            # try:
            bfs = main(shot, plunge, showplot=False, min_npeaks=5)
            # except Exception as e:
            #    print("Failed")
            #    print(e)
            #    continue
            min_r = max(bfs.rcp_r.min(), bfs.blob_r.min())
            max_r = min(bfs.rcp_r.max(), bfs.blob_r.max())
            blob_mask = np.logical_and(bfs.blob_r > min_r, bfs.blob_r < max_r)
            if type(smooth_window) == type(None):
                # f_vr = interp1d(bfs.rcp_r, bfs.vr_solved)
                # f_vr = interp1d(bfs.rcp_r, bfs.vr_solved2)
                f_vr = interp1d(bfs.rcp_r, bfs.vr_solved3)
                f_lambda_ne = interp1d(bfs.rcp_r, bfs.lambda_ne)
            else:
                # f_vr = interp1d(bfs.rcp_r, savgol_filter(bfs.vr_solved, smooth_window, 2))
                # f_vr = interp1d(bfs.rcp_r, savgol_filter(bfs.vr_solved2, smooth_window, 2))
                f_vr = interp1d(bfs.rcp_r, savgol_filter(bfs.vr_solved3, smooth_window, 2))
                f_lambda_ne = interp1d(bfs.rcp_r, savgol_filter(bfs.lambda_ne, smooth_window, 2))
            f_conn = interp1d(bfs.conns_r, bfs.conns_l_otf + bfs.conns_l_itf)
            f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
            f_E = interp1d(bfs.rcp_r, bfs.E_solved)
            f_M = interp1d(bfs.rcp_r, bfs.rcp_m)
            if plot_data == "vr":
                x = np.append(x, f_vr(bfs.blob_r[blob_mask]))
                y = np.append(y, bfs.blob_vr_time_avg[blob_mask])
                y_ExB = np.append(y_ExB, f_E(bfs.blob_r[blob_mask]) / 1.5)
            elif plot_data == "lambda_ne":
                x = np.append(x, bfs.solved_lambda_ne[blob_mask])
                y = np.append(y, f_lambda_ne(bfs.blob_r[blob_mask]))
            # y = np.append(y, bfs.blob_vr[blob_mask])
            r = np.append(r, bfs.blob_rmrs[blob_mask])
            # r = np.append(r, bfs.blob_r[blob_mask])
            te = np.append(te, bfs.blob_te[blob_mask])
            ne = np.append(ne, bfs.blob_ne[blob_mask])
            drad = np.append(drad, bfs.blob_drad[blob_mask])
            conn = np.append(conn, f_conn(bfs.blob_r[blob_mask]))
            lambda_ne = np.append(lambda_ne, f_lambda_ne(bfs.blob_r[blob_mask]))
            dvdr = np.append(dvdr, np.full(np.sum(blob_mask), bfs.dvdr))
            M = np.append(M, f_M(bfs.blob_r[blob_mask]))
            vr = np.append(vr, bfs.blob_vr[blob_mask])
            f = np.append(f, bfs.blob_f[blob_mask])
            nesep = np.append(nesep, np.full(np.sum(blob_mask), bfs.nesep))
            twidth = np.append(twidth, bfs.blob_twidth[blob_mask])
            f_isat = interp1d(bfs.rcp_r, bfs.rcp_isat, bounds_error=False)
            isat = np.append(isat, f_isat(bfs.blob_r[blob_mask]))

            if plot_opt == 1:
                if count < 10:
                    ax1.scatter(f_vr(bfs.blob_r[blob_mask]), bfs.blob_vr_time_avg[blob_mask],
                                label="{}_{}".format(shot, plunge))
                else:
                    ax1.scatter(f_vr(bfs.blob_r[blob_mask]), bfs.blob_vr_time_avg[blob_mask],
                                label="{}_{}".format(shot, plunge), marker="^")

            # For this plot we probably want to apply the connection mask as well.
            elif plot_opt == 2:

                blob_mask = np.logical_and(blob_mask, f_conn(bfs.blob_r) > 0)
                blob_mask = np.logical_and(blob_mask, bfs.blob_rmrs > 3)
                blob_mask = np.logical_and(blob_mask, f_isat(bfs.blob_r) > 0.025)

                if plot_data == "vr":
                    opt2x = np.append(opt2x, f_vr(bfs.blob_r[blob_mask].mean()))
                    opt2y = np.append(opt2y, bfs.blob_vr_time_avg[blob_mask].mean())
                    opt2c = np.append(opt2c, f_lambda_ne(bfs.blob_r[blob_mask]).mean())
                elif plot_data == "lambda_ne":
                    # opt2y = np.append(opt2y, f_lambda_ne(bfs.blob_r[blob_mask].mean()))
                    # opt2x = np.append(opt2x, bfs.solved_lambda_ne2[blob_mask].mean())
                    # opt2c = np.append(opt2c, bfs.neut_dens2[blob_mask].mean())
                    opt2y = np.append(opt2y, f_lambda_ne(bfs.blob_r[blob_mask]))
                    opt2x = np.append(opt2x, bfs.solved_lambda_ne2[blob_mask])
                    # opt2c = np.append(opt2c, bfs.neut_dens2[blob_mask])
                    # opt2c = np.append(opt2c, 235 - bfs.blob_r[blob_mask])  # Distance from wall/poloidal bumper.
                    # opt2c = np.append(opt2c, bfs.blob_ne[blob_mask])
                    # opt2c = np.append(opt2c, 2.75e-14 * bfs.blob_ne[blob_mask] * f_conn(bfs.blob_r[blob_mask]) /
                    #                   np.square(bfs.blob_te[blob_mask]))
                    # opt2c = np.append(opt2c, bfs.blob_rmrs[blob_mask])

                    # The ei-collision frequency calculation. From page 197 in Friedberg and Eq. 20 in Myra PoP 2006.
                    # f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
                    vi = np.sqrt(3 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    vte = np.sqrt(2 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    nu_ei = (1 / (4 * np.pi) * bfs.blob_ne[blob_mask] * elec ** 4 / (eps ** 2 * me * md) * ln_lambda) * \
                            (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))  # 1/s
                    elec_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / me
                    ion_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / md
                    ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * bfs.blob_te[blob_mask]) / (
                            elec * bfs.blob_b[blob_mask])
                    norm_coll = nu_ei * f_conn_otf(bfs.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)
                    opt2c = np.append(opt2c, norm_coll)

                    if shot in [184527, 187111]:
                        opt2m = np.append(opt2m, np.full(sum(blob_mask), "*"))
                    else:
                        opt2m = np.append(opt2m, np.full(sum(blob_mask), "^"))

                elif plot_data == "lambda_ne_neut":

                    opt2y = np.append(opt2y, f_lambda_ne(bfs.blob_r[blob_mask]))
                    opt2x = np.append(opt2x, bfs.solved_lambda_ne4[blob_mask])
                    opt2c = np.append(opt2c, bfs.blob_rmrs[blob_mask])

                    if shot in [184527, 187111]:
                        opt2m = np.append(opt2m, np.full(sum(blob_mask), "*"))
                    else:
                        opt2m = np.append(opt2m, np.full(sum(blob_mask), "^"))

                # Interested in plotting the needed neutral density against the difference in the predicted lambda_ne
                # and measured.
                elif plot_data == "neut_dens":
                    diff = f_lambda_ne(bfs.blob_r[blob_mask]) / bfs.solved_lambda_ne2[blob_mask]
                    # opt2y = np.append(opt2y, bfs.neut_dens2[blob_mask].mean())
                    # opt2x = np.append(opt2x, diff.mean())
                    # opt2c = np.append(opt2c, bfs.blob_vr_time_avg[blob_mask].mean())
                    opt2y = np.append(opt2y, bfs.neut_dens2[blob_mask])
                    opt2x = np.append(opt2x, diff)
                    # opt2c = np.append(opt2c, bfs.blob_vr_time_avg[blob_mask])
                    opt2c = np.append(opt2c, bfs.blob_ne[blob_mask])

                # Slightly different plot Where we are plotting the radial blob velocity against the neutral density,
                # and colorcoding by the lambda. Could be presented as a contour plot best.
                elif plot_data == "category":
                    # if shot in [190484, 190485, 190486]:
                    opt2x = np.append(opt2x, bfs.neut_dens2[blob_mask])
                    opt2y = np.append(opt2y, bfs.blob_vr_time_avg[blob_mask])
                    opt2c = np.append(opt2c, f_lambda_ne(bfs.blob_r[blob_mask]) / f_conn(bfs.blob_r[blob_mask]))

                elif plot_data == "myra":

                    # Need to calculate the normalized collisionality and blob size. First calculate collisionality.
                    # f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
                    vi = np.sqrt(3 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    vte = np.sqrt(2 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    nu_ei = (1 / (4 * np.pi) * bfs.blob_ne[blob_mask] * elec ** 4 / (eps ** 2 * me * md) * ln_lambda) * \
                            (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))  # 1/s
                    elec_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / me
                    ion_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / md
                    ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * bfs.blob_te[blob_mask]) / (
                            elec * bfs.blob_b[blob_mask])
                    norm_coll = nu_ei * f_conn_otf(bfs.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)

                    # Now calculate the normalized radius.
                    ab = bfs.blob_drad[blob_mask] / 2 / 100  # cm to m
                    R = bfs.blob_r[blob_mask] / 100  # cm to m
                    L = f_conn_otf(bfs.blob_r[blob_mask])
                    ahat = ab * np.power(R, 1 / 5) / (np.power(L, 2 / 5) * np.power(ion_gyrorad, 4 / 5))

                    norm_rad = np.power(ahat, 5 / 2)

                    # Now the normalized velocity.
                    vr = bfs.blob_vr[blob_mask]
                    cs = bfs.blob_cs[blob_mask]
                    norm_vr = vr / (np.power(2 * L * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs)

                    # Append to plot data.
                    ahats = np.append(ahats, ahat)
                    opt2x = np.append(opt2x, norm_rad)
                    opt2y = np.append(opt2y, norm_coll)
                    opt2c = np.append(opt2c, norm_vr)

                elif plot_data == "instability":

                    # Need to calculate the normalized collisionality and blob size. First calculate collisionality.
                    # f_conn_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
                    vi = np.sqrt(3 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    vte = np.sqrt(2 * elec * bfs.blob_te[blob_mask] / md)  # m/s
                    nu_ei = (1 / (4 * np.pi) * bfs.blob_ne[blob_mask] * elec ** 4 / (eps ** 2 * me * md) * ln_lambda) * \
                            (1 / (np.power(vi, 3) + 1.3 * np.power(vte, 3)))  # 1/s
                    elec_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / me
                    ion_gyrofreq = 1.609e-19 * bfs.blob_b[blob_mask] / md
                    ion_gyrorad = np.sqrt(2 * 2.014 * 1.66e-27 * elec * bfs.blob_te[blob_mask]) / (
                            elec * bfs.blob_b[blob_mask])
                    norm_coll = nu_ei * f_conn_otf(bfs.blob_r[blob_mask]) / (elec_gyrofreq * ion_gyrorad)

                    # Now calculate the normalized radius.
                    ab = bfs.blob_drad[blob_mask] / 2 / 100  # cm to m
                    R = bfs.blob_r[blob_mask] / 100  # cm to m
                    L = f_conn_otf(bfs.blob_r[blob_mask])
                    ahat = ab * np.power(R, 1 / 5) / (np.power(L, 2 / 5) * np.power(ion_gyrorad, 4 / 5))

                    norm_rad = np.power(ahat, 5 / 2)

                    # Now the normalized velocity.
                    vr = bfs.blob_vr[blob_mask]
                    cs = bfs.blob_cs[blob_mask]
                    norm_vr = vr / (np.power(2 * L * np.square(ion_gyrorad) / np.power(R, 3), 1 / 5) * cs)

                    # Calculate the blob density above the background value.
                    f_rcp_ne = interp1d(bfs.rcp_r, bfs.rcp_ne)
                    deltn_n = (bfs.blob_ne[blob_mask] - f_rcp_ne(bfs.blob_r[blob_mask])) / f_rcp_ne(
                        bfs.blob_r[blob_mask])

                    # Calculate the ion-neutral collision (ionization) rate. Note: This may not be correct since we
                    # want *ion*-neutral, and ionization is an electron-neutral thing.
                    nu_iz = bfs.neut_dens2[blob_mask] * bfs.blob_ioniz_rate_coeffs[blob_mask]

                    # Calculate the ion-neutral collision (CX) rate.
                    nu_cx = bfs.neut_dens2[blob_mask] * bfs.blob_cx_rate_coeffs[blob_mask]

                    # Now calculate the secondary instability growth rate, Eq. 5.1.13 from Theiler's thesis.
                    B = bfs.blob_b[blob_mask]
                    larmor = np.sqrt(bfs.blob_te[blob_mask] * md * elec) / (elec * B)
                    yint = np.sqrt(2 / (R * ab)) * cs
                    instab_ratio = np.sqrt(2 * ab / R) * cs * deltn_n / vr - (1 / (np.square(larmor) * L)) \
                                   * np.sqrt(R / 2) * np.power(ab, 5 / 2) - nu_cx * np.sqrt(R * ab) / (np.sqrt(2) * cs)
                    yinstab = instab_ratio * yint

                    window = int(len(instab_ratio) / 2)
                    if window % 2 == 0:
                        window += 1
                    if window <= 2:
                        opt2x = np.append(opt2x, instab_ratio)
                        # opt2x = np.append(opt2x, yinstab)
                        opt2y = np.append(opt2y, norm_vr)
                        opt2c = np.append(opt2c, yinstab)
                        # opt2c = np.append(opt2c, norm_rad)
                    else:
                        opt2x = np.append(opt2x, savgol_filter(instab_ratio, window, 2))
                        # opt2x = np.append(opt2x, savgol_filter(yinstab, window, 2))
                        opt2y = np.append(opt2y, savgol_filter(norm_vr, window, 2))
                        opt2c = np.append(opt2c, savgol_filter(yinstab, window, 2))
                        # opt2c = np.append(opt2c, savgol_filter(norm_rad, window, 2))

                # opt2c = np.append(opt2c, 4.2e13 * np.power(bfs.blob_te[blob_mask].mean(), 3/2) / bfs.blob_ne[blob_mask].mean())
                labels = np.append(labels, np.full(np.sum(blob_mask), "{}_{}".format(shot, plunge)))
                # if count < 10:
                #     ax1.scatter(f_vr(bfs.blob_r[blob_mask]).mean(), bfs.blob_vr_time_avg[blob_mask].mean(),
                #         label="{}_{}".format(shot, plunge))
                # else:
                #     ax1.scatter(f_vr(bfs.blob_r[blob_mask]).mean(), bfs.blob_vr_time_avg[blob_mask].mean(),
                #         label="{}_{}".format(shot, plunge), marker="^")

            # Plots of the plasma density normalized to the separatrix value.
            # nenorm = bfs.rcp_ne / bfs.nesep
            # if count < 10:
            #     ax2.plot(bfs.rcp_rho, nenorm, label="{}_{}".format(shot, plunge))
            # else:
            #     ax2.plot(bfs.rcp_rho, nenorm, label="{}_{}".format(shot, plunge), linestyle="--")

            count += 1

    # Create colorscale according to measurement. Point is to see if the trend breaks down at certain values.
    sol_coll = 1e-16 * ne * conn / np.square(te)
    conn_mask = conn > 20

    from matplotlib.colors import LogNorm
    if plot_opt == 3:
        ax1.plot([0, 500], [0, 500], color="k", linestyle="--")
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=sol_coll[conn_mask], cmap="magma", norm=LogNorm(vmin=50, vmax=1000))
        # sc = ax1.scatter(x[conn_mask], vr[conn_mask], c=sol_coll[conn_mask], cmap="magma", norm=LogNorm(vmin=50, vmax=1000))
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=r[conn_mask], cmap="magma")
        # sc = ax1.scatter(x, y, c=lambda_ne, cmap="magma", vmin=0, vmax=10)
        # sc = ax1.scatter(x, y, c=dvdr, cmap="magma", vmin=-700, vmax=0)
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=drad[conn_mask], cmap="magma", vmin=0, vmax=4)
        sc = ax1.scatter(x[conn_mask], y[conn_mask], c=conn[conn_mask], cmap="magma")
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=te[conn_mask], cmap="magma")
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=M[conn_mask], cmap="magma", vmin=-0.2, vmax=0.2)
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=f[conn_mask], cmap="magma", vmin=1000, vmax=3000)
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=nesep[conn_mask], cmap="magma")
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=twidth[conn_mask], cmap="magma", norm=LogNorm())
        # sc = ax1.scatter(x[conn_mask], y[conn_mask], c=isat[conn_mask], cmap="magma")

        cbar = fig.colorbar(sc, ax=ax1)
        # cbar.set_label("SOL Collisionality")
        # cbar.set_label("R-Rsep (cm)")
        cbar.set_label("Connection Length (m)")
        # cbar.set_label("lambda_ne (cm)")
        # cbar.set_label("dvr / dr")
        # cbar.set_label("Drad (cm)")
        # cbar.set_label("Blob Frequency (Hz)")
        # cbar.set_label("nesep (m-3)")
        # cbar.set_label("Tblob (s)")
        # cbar.set_label("Isat (A)")

        if plot_data == "vr":
            ax1.set_xlabel("Simple SOL vr", fontsize=14)
            ax1.set_ylabel("Time-averaged blob vr", fontsize=14)
            ax1.set_xlim([0, 100])
            ax1.set_ylim([0, 60])
        elif plot_data == "lambda_ne":
            ax1.set_xlabel(r"Predicted $\mathdefault{\lambda_{ne}}$ (cm)", fontsize=16)
            ax1.set_ylabel(r"Measured $\mathdefault{\lambda_{ne}}$ (cm)", fontsize=16)
            ax1.tick_params(axis="both", which="major", labelsize=14)
            ax1.set_xlim([0, 20])
            ax1.set_ylim([0, 20])

    # Taken from https://github.com/matplotlib/matplotlib/issues/11155 to allow individual markers per point.
    def mscatter(x, y, ax=None, m=None, **kw):
        import matplotlib.markers as mmarkers
        if not ax: ax = plt.gca()
        sc = ax.scatter(x, y, **kw)
        if (m is not None) and (len(m) == len(x)):
            paths = []
            for marker in m:
                if isinstance(marker, mmarkers.MarkerStyle):
                    marker_obj = marker
                else:
                    marker_obj = mmarkers.MarkerStyle(marker)
                path = marker_obj.get_path().transformed(
                    marker_obj.get_transform())
                paths.append(path)
            sc.set_paths(paths)
        return sc

    if plot_opt == 2:
        if plot_data == "vr":
            ax1.plot([0, 500], [0, 500], color="k", linestyle="--")
            sc = ax1.scatter(opt2x, opt2y, c=opt2c, cmap="inferno", vmin=0, vmax=15, zorder=15, marker="^",
                             edgecolors="k", s=75)
            cbar = fig.colorbar(sc, ax=ax1)
            cbar.set_label(r"$\mathdefault{\lambda_{ne}}$ (cm)", fontsize=14)
            ax1.set_xlabel("Simple SOL vr", fontsize=14)
            ax1.set_ylabel("Time-averaged blob vr", fontsize=14)
            ax1.set_xlim([0, 100])
            ax1.set_ylim([0, 60])
            ax1.grid(zorder=5)
        elif plot_data == "lambda_ne":
            ax1.plot([0, 500], [0, 500], color="k", linestyle="--")
            sc = mscatter(opt2x, opt2y, ax=ax1, m=opt2m, c=opt2c, cmap="inferno", zorder=15, marker="^", edgecolors="k",
                          s=75)
            # sc = ax1.scatter(opt2x, opt2y, c=opt2c, cmap="inferno", vmin=2, vmax=10, zorder=15, marker="^",
            #                  edgecolors="k", s=75)
            cbar = fig.colorbar(sc, ax=ax1)
            # cbar.set_label(r"Neutral Density ($\mathdefault{m^{-3}}$)", fontsize=14)
            cbar.set_label(r"$\mathdefault{R-R_{sep}\ (cm)}$", fontsize=16)
            ax1.set_xlabel(r"Predicted $\mathdefault{\lambda_{ne}}$ (cm)", fontsize=16)
            ax1.set_ylabel(r"Measured $\mathdefault{\lambda_{ne}}$ (cm)", fontsize=16)
            ax1.tick_params(axis="both", which="major", labelsize=14)
            ax1.set_xlim([0, 8])
            ax1.set_ylim([0, 14])
            ax1.grid(zorder=5)

            # A separate plot of the Z data against the ratio of the lambdas.
            fig2, ax2 = plt.subplots(1, 1, figsize=(5, 4))
            ax2.scatter(opt2c, opt2y / opt2x, zorder=15, marker="^", edgecolors="k", s=75)
            ax2.grid(zorder=5)
            ax2.set_xlabel(r"$\mathdefault{R-R_{wall}\ (cm)}$", fontsize=14)
            ax2.set_ylabel(r"$\mathdefault{\lambda_{ne}^{meas} / \lambda_{ne}^{pred}}$", fontsize=14)
            fig2.tight_layout()
            fig2.show()

        elif plot_data == "lambda_ne_neut":
            ax1.plot([0, 500], [0, 500], color="k", linestyle="--")
            sc = mscatter(opt2x, opt2y, ax=ax1, m=opt2m, c=opt2c, cmap="inferno", zorder=15, marker="^", edgecolors="k",
                          s=75)
            # sc = ax1.scatter(opt2x, opt2y, c=opt2c, cmap="inferno", vmin=2, vmax=10, zorder=15, marker="^",
            #                  edgecolors="k", s=75)
            cbar = fig.colorbar(sc, ax=ax1)
            # cbar.set_label(r"Neutral Density ($\mathdefault{m^{-3}}$)", fontsize=14)
            cbar.set_label(r"$\mathdefault{R-R_{sep}\ (cm)}$", fontsize=16)
            ax1.set_xlabel(r"Predicted $\mathdefault{\lambda_{ne}}$ incl. LLAMA (cm)", fontsize=16)
            ax1.set_ylabel(r"Measured $\mathdefault{\lambda_{ne}}$ (cm)", fontsize=16)
            ax1.tick_params(axis="both", which="major", labelsize=14)
            ax1.set_xlim([0, 8])
            ax1.set_ylim([0, 14])
            ax1.grid(zorder=5)

            # A separate plot of the Z data against the ratio of the lambdas.
            fig2, ax2 = plt.subplots(1, 1, figsize=(5, 4))
            ax2.scatter(opt2c, opt2y / opt2x, zorder=15, marker="^", edgecolors="k", s=75)
            ax2.grid(zorder=5)
            ax2.set_xlabel(r"$\mathdefault{R-R_{wall}\ (cm)}$", fontsize=14)
            ax2.set_ylabel(r"$\mathdefault{\lambda_{ne}^{meas} / \lambda_{ne}^{pred}}$", fontsize=14)
            fig2.tight_layout()
            fig2.show()

        elif plot_data == "neut_dens":
            ax1.plot([0, 500], [0, 500], color="k", linestyle="--")
            # sc = ax1.scatter(opt2x, opt2y, c=opt2c, cmap="inferno", vmin=0, vmax=50, zorder=15, marker="^",
            #                  edgecolors="k", s=75)
            sc = ax1.scatter(opt2x, opt2y, c="tab:red", zorder=15, marker="^",
                             edgecolors="k", s=75)
            # cbar = fig.colorbar(sc, ax=ax1)
            # cbar.set_label(r"$\mathdefault{v_{r}}$ (m/s)", fontsize=14)
            ax1.set_xlabel(r"$\mathdefault{\lambda_{ne}^{meas} / \lambda_{ne}^{pred}}$", fontsize=16)
            ax1.set_ylabel(r"$\mathdefault{n_n\ (m^{-3})}$", fontsize=16)
            ax1.set_xlim([0, 10])
            ax1.set_ylim([5e15, 1e18])
            ax1.set_yscale("log")
            ax1.grid(zorder=5, which="both", alpha=0.4)
            ax1.tick_params(axis="both", which="major", labelsize=14)

        elif plot_data == "category":
            mask = ~np.isnan(opt2x)
            print(opt2x[mask])
            print(opt2y[mask])
            print(opt2c[mask])
            # ax1.tricontourf(opt2x[mask], opt2y[mask], opt2c[mask], cmap="inferno")
            sc = ax1.scatter(opt2x, opt2y, c=opt2c, cmap="magma", zorder=15, vmin=0, vmax=0.2)
            cbar = fig.colorbar(sc, ax=ax1)
            cbar.set_label("lambda / conn")
            ax1.set_xlabel(r"$\mathdefault{n_n\ (m^{-3})}$", fontsize=16)
            ax1.set_ylabel("Time-averaged blob vr", fontsize=14)
            ax1.set_xscale("log")
            ax1.grid(zorder=5)

        elif plot_data == "myra":
            eps = 0.3  # Rough approximation of the flux expansion factor.
            ax1.plot([0, 50], [0, 50], color="k")
            ax1.plot([0, 1 / eps], [0, 1], color="k")
            ax1.plot([1 / eps, 50], [1, 1], color="k")
            ax1.plot([1 / eps, 1 / eps], [0, 1], color="k")
            sc = ax1.scatter(opt2x, opt2y, c=opt2c, marker="^", cmap="magma", zorder=15, edgecolors="k", s=75)
            cbar = fig.colorbar(sc, ax=ax1)
            ax1.grid(zorder=5, which="both", alpha=0.3)
            cbar.set_label(r"$\mathdefault{\hat{v_r}}$")
            ax1.set_xlabel(r"$\mathdefault{\Theta}$", fontsize=16)
            ax1.set_ylabel(r"$\mathdefault{\Lambda_{OMP}}$", fontsize=16)
            ax1.set_yscale("log")
            ax1.set_xscale("log")
            ax1.set_xlim([0.1, 50])
            ax1.set_ylim([0.1, 50])

        elif plot_data == "instability":
            ax1.axvline(0.0, color="k", linestyle="--")
            sc = ax1.scatter(opt2x, opt2y, c=opt2c, marker="^", cmap="coolwarm", zorder=15, edgecolors="k", s=75,
                             vmin=-0.3e7, vmax=0.3e7)
            # sc = ax1.scatter(opt2x, opt2y, c=opt2c, marker="^", cmap="magma", zorder=15, edgecolors="k", s=75,
            #                  vmin=0, vmax=3)
            # sc = ax1.scatter(opt2x, opt2c, marker="^", color="tab:red", zorder=15, edgecolors="k", s=75)
            # cbar = fig.colorbar(sc, ax=ax1)
            ax1.grid(zorder=5, which="both", alpha=0.3)
            cbar.set_label(r"$\gamma_{inst}$ (s)", fontsize=14)
            # cbar.set_label(r"$\mathdefault{\hat{a}}$", fontsize=14)
            ax1.set_xlabel(r"$\gamma_{inst}$ \ $\gamma_{int}$", fontsize=16)
            # ax1.set_xlabel(r"$\gamma_{inst}$ (s)", fontsize=16)
            ax1.set_ylabel(r"$\mathdefault{\hat{v_r}}$", fontsize=16)
            # ax1.set_ylabel(r"$\mathdefault{\hat{a}}$", fontsize=16)
            # ax1.set_ylabel(r"$\mathdefault{\Theta}$", fontsize=16)
            # ax1.set_xlim([-22, 7])
            # ax1.set_xlim([-8e6, 3e6])
            # ax1.set_ylim([0, 0.8])
            # ax1.set_ylim([0, 3.5])

        for i in range(0, len(opt2x)):
            ax1.annotate(labels[i], (opt2x[i], opt2y[i]))

    # ax1.legend()

    # ax2.scatter(drad[conn_mask], sol_coll[conn_mask], c=y[conn_mask], cmap="magma")
    # ax2.set_xlabel("Blob width (cm)")
    # ax2.set_ylabel("SOL collisionality")
    # ax2.grid()
    # #ax2.set_xscale("log")
    # ax2.legend()
    # ax2.set_ylabel("ne/nesep")
    # ax2.set_xlabel("Rho")

    fig.tight_layout()
    fig.show()

    # Additional plot of the normalized velocity against the ahat, similar to Cedric's PoP 2018 paper.
    if plot_data == "myra":
        fig, ax1 = plt.subplots(figsize=(5, 4))
        ax1.plot(np.linspace(0, 5, 50), np.sqrt(np.linspace(0, 5, 50)), color="k", linestyle="--", lw=3)
        ax1.scatter(ahats, opt2c, marker="^", s=75, color="tab:red", edgecolors="k")
        ax1.set_ylabel(r"$\mathdefault{\hat{v_r}}$", fontsize=16)
        ax1.set_xlabel(r"$\mathdefault{\hat{a}}$", fontsize=16)
        ax1.set_ylim([0, 2])
        ax1.set_xlim([0, 5])
        fig.tight_layout()
        fig.show()

        fig, ax1 = plt.subplots(figsize=(5, 4))
        # ax1.plot(np.linspace(0, 5, 50), np.sqrt(np.linspace(0, 5, 50)), color="k", linestyle="--", lw=3)
        x = opt2y / np.square(ahats)
        ax1.plot([0, 40], [0, 40])
        ax1.scatter(x, opt2c, marker="^", s=75, color="tab:red", edgecolors="k")
        ax1.set_ylabel(r"$\mathdefault{\hat{v_r}}$", fontsize=16)
        ax1.set_xlabel(r"$\mathdefault{\Lambda_{OMP}/\hat{a}^2}$", fontsize=16)
        ax1.set_ylim([0, 1])
        ax1.set_xlim([0, 40])
        fig.tight_layout()
        fig.show()


def fit_lambda_ne(shot, plunge, minidx=0, maxidx=9999, mach_mult=1):
    # Use the already coded path to save time.
    bfs = main(shot, plunge, showplot=False)

    # Easier to do a linear fit of ln(ne). Pull that out and then create a subset of data that will be fit.
    r = bfs.rcp_r
    lnne = np.log(bfs.rcp_ne)
    r_fit = r[minidx:maxidx]
    lnne_fit = lnne[minidx:maxidx]
    cs_fit = bfs.rcp_cs[minidx:maxidx]
    te_fit = bfs.rcp_te[minidx:maxidx]
    ne_fit = bfs.rcp_ne[minidx:maxidx]

    # Perform fit.
    try:
        z = np.polyfit(r_fit, lnne_fit, 1)
    except:
        print("Error in linear fit!")
        print(" r_fit: {}".format(r_fit))
        print(" lnne_fit: {}".format(lnne_fit))
    p = np.poly1d(z)
    lnne_fitline = p(r_fit)
    lambda_ne = -1 / z[0]
    print("lambda_ne = {:.2f} cm".format(lambda_ne))
    print("Fit range: {:.2f}-{:.2f} cm".format(r_fit.min(), r_fit.max()))

    # Get weighted average of blob traits in this fit range.
    blob_mask = np.logical_and(bfs.blob_r > r_fit.min(), bfs.blob_r < r_fit.max())
    blob_vr_fit = bfs.blob_vr_time_avg[blob_mask]
    blob_npeaks_fit = bfs.blob_npeaks[blob_mask]
    avg_vr = np.average(blob_vr_fit, weights=blob_npeaks_fit)
    print("Avg. vr = {:.2f} m/s".format(avg_vr))

    # Copied from above, estimate a C value.
    f_itf = interp1d(bfs.conns_r, bfs.conns_l_itf)
    f_otf = interp1d(bfs.conns_r, bfs.conns_l_otf)
    rcp_lconn_itf = f_itf(r_fit)
    rcp_lconn_otf = f_otf(r_fit)
    rcp_s = rcp_lconn_itf
    rcp_smax = rcp_lconn_otf + rcp_lconn_itf

    # Again, copied from above. This is an equally valid L definition to use in C.
    meas_machs = bfs.rcp_m[minidx:maxidx] * mach_mult
    L = np.zeros(len(r_fit))
    for i in range(0, len(r_fit)):
        meas_mach = meas_machs[i]
        if meas_mach > 0:
            slope = (meas_mach - 1) / (rcp_s[i] - rcp_smax[i])
        else:
            slope = (meas_mach + 1) / (rcp_s[i] - 0)
        stag_s = -meas_mach / slope + rcp_s[i]
        it_conn = stag_s
        ot_conn = rcp_smax[i] - stag_s
        if meas_mach > 0:
            L[i] = ot_conn
        else:
            L[i] = it_conn
        s = np.abs(rcp_s[i] - stag_s)

    C1 = (np.pi / 2 - 1) * cs_fit / (rcp_smax / 2)
    print("[1] Avg. C = {:.2f}".format(C1.mean()))
    C2 = (np.pi / 2 - 1) * cs_fit / L
    print("[2] Avg. C = {:.2f}".format(C2.mean()))
    # print("    C1 --> C2")
    # for i in range(0, len(L)):
    #    print(" {:.2f} --> {:.2f}".format(rcp_smax[i]/2, L[i]))

    # The average neutral density estimate.
    oa = openadas.OpenADAS()
    rate_df = oa.read_rate_coef_unres("/Users/zamperini/My Drive/Research/Data/openadas/scd96_h.dat")
    ioniz_rate_coef = oa.get_rate_coef(rate_df, te_fit.mean(), ne_fit.mean())
    neut_dens = ((np.pi / 2 - 1) * cs_fit.mean() / (rcp_smax.mean() / 2) - avg_vr /
                 (lambda_ne / 100)) / ioniz_rate_coef
    print("Average neutral density: {:.2e}".format(neut_dens))
    print("Average ne: {:.2e}".format(ne_fit.mean()))

    # Again but using the newer equations.
    lambdas = np.zeros(len(r_fit))
    neuts = np.zeros(len(r_fit))
    for i in range(0, len(r_fit)):
        if rcp_s[i] > rcp_smax[i] / 2:
            s = rcp_s[i] - rcp_smax[i] / 2
        else:
            s = rcp_smax[i] / 2 - rcp_s[i]
        L = rcp_smax[i] / 2

        # To calculate n0 for the parallel, take a given set of (s, ne) and using the linear Mach assumption just
        # back out what n0 is.
        tmp_M = rcp_s[i] / rcp_smax[i] - 1
        n0 = ne_fit[i] * (1 + tmp_M ** 2)

        dnds = n0 * 2 * L ** 2 * (L - s) / (2 * L ** 2 - 2 * L * s + s ** 2) ** 2
        vpar = (s / L - 1) * cs_fit[i]
        dvds = cs_fit[i] / L
        dvdr = 0
        tmp = dnds * vpar / (ne_fit[i] * blob_vr_fit.mean()) + 1 / blob_vr_fit.mean() * dvds + 1 / \
              blob_vr_fit.mean() * dvdr
        lambdas[i] = 1 / tmp

        # Again calculate this.
        # print("term 1 = {}".format(vpar / ne_fit[i] * dnds))
        # print("term 2 = {}".format(dvds))
        # print("term 3 = {}".format(dvdr))
        # print("term 4 = {}".format(blob_vr_fit.mean() / (lambda_ne / 100)))
        tmp = vpar / ne_fit[i] * dnds + dvds + dvdr - blob_vr_fit.mean() / (lambda_ne / 100)
        ioniz_rate_coef = oa.get_rate_coef(rate_df, te_fit[i], ne_fit[i])
        neuts[i] = tmp / ioniz_rate_coef

    print("From newer derivation")
    print(" Average lambda_ne: {:.2} cm".format(lambdas.mean() * 100))
    print(" Average neutral density: {:.2e}".format(neuts.mean()))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    ax1.scatter(r, lnne, color="k", s=35)
    ax1.scatter(r_fit, lnne_fit, color="r", s=20)
    ax1.plot(r_fit, lnne_fitline, color="r")
    ax1.grid()
    ax1.set_xlabel("R (cm)")
    ax1.set_ylabel("ln(ne)")
    ax1.set_title("{}_{}".format(shot, plunge))
    ax2.axvline(r_fit.min(), color="k")
    ax2.axvline(r_fit.max(), color="k")
    ax2.plot(bfs.conns_r, bfs.conns_l_itf, label="ITF")
    ax2.plot(bfs.conns_r, bfs.conns_l_otf, label="OTF")
    ax2.legend()
    ax2.set_xlabel("R (cm)")
    ax2.set_ylabel("Connection Length (m)")
    ax2.set_yscale("log")
    fig.tight_layout()
    fig.show()

    return bfs


if __name__ == "__main__":
    if sys.argv[1] == "all":
        if len(sys.argv) == 3:
            main_all(int(sys.argv[2]))
        else:
            main_all()
    elif sys.argv[1] == "fit":
        if len(sys.argv) != 6:
            print("Incorrect usage! Example usage:")
            print("  python BlobbyFarSOL.py fit shot plunge minidx maxidx")
            print("  python BlobbyFarSOL.py fit 190440 1 2 18")
        shot = int(sys.argv[2])
        plunge = int(sys.argv[3])
        minidx = int(sys.argv[4])
        maxidx = int(sys.argv[5])
        if shot in [184527, 187111, 167193, 167195]:
            mach_mult = -1.0
        else:
            mach_mult = 1.0
        bfs = fit_lambda_ne(shot, plunge, minidx, maxidx, mach_mult)

    else:
        shot = int(sys.argv[1])
        plunge = int(sys.argv[2])
        bfs = main(shot, plunge)
