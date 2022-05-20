# This script compares Zeff for a set of DIVIMP runs considering graphite
# wall and another considering SiC walls.
import oedge_plots
import numpy as np


grap_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-grap-shelf.nc"
sic_c_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-sic-c-shelf.nc"
sic_si_path = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-mm-sic-si-shelf.nc"

print("Loading runs...")
grap = oedge_plots.OedgePlots(grap_path)
sic_c = oedge_plots.OedgePlots(sic_c_path)
sic_si = oedge_plots.OedgePlots(sic_si_path)


def zeff_grap(op):

    # To calculate Zeff, it will be the weighted average of nz_charge / ne,
    # where the weights are the proportion of each charge state, i.e.,
    # ddlims_charge / ddlims_sum.
    ne = op.read_data_2d("KNBS")
    #ddlims_sum = op.read_data_2d("DDLIMS", charge="all", scaling=op.absfac)
    cion = int(op.nc.variables["CION"][:])

    #tot_dens = ddlims_sum + ne
    zeff = ne.copy()  # Approximating nD ~ ne here

    # Total density of ions is sum of each charge state plus ne.
    #tot_dens = ne.copy()

    #ddlims_charges = []
    for charge in range(0, cion):

        # charge + 1 so that we start at C1+ and not the neutrals (C0+).
        impdens = op.read_data_2d("DDLIMS", charge=charge+1, scaling=op.absfac)
        #ddlims_charges.append(impdens)
        #tot_dens += impdens
        zeff += impdens * (charge + 1)**2

    #for i in range(0, len(ne)):
    #    for charge in range(0, cion):
    #        zeff[i] += ddlims_charges[charge][i] * (charge + 1)

    # Right now zeff is the total charge/m-3. Divide by tot_dens to get zeff.
    #zeff = zeff / tot_dens
    return zeff / ne

def zeff_sic(op_c, op_si):

    # To calculate Zeff, it will be the weighted average of nz_charge / ne,
    # where the weights are the proportion of each charge state, i.e.,
    # ddlims_charge / ddlims_sum.
    ne = op_c.read_data_2d("KNBS")
    #ddlims_sum_c = op_c.read_data_2d("DDLIMS", charge="all", scaling=op_c.absfac)
    #ddlims_sum_si = op_si.read_data_2d("DDLIMS", charge="all", scaling=op_si.absfac)
    cion_c = int(op_c.nc.variables["CION"][:])
    cion_si = int(op_si.nc.variables["CION"][:])

    #tot_dens = ddlims_sum_si + ddlims_sum_c + ne
    zeff = ne.copy()  # Approximating nD ~ ne here

    #tot_dens = ne.copy()

    #ddlims_charges_c = []
    for charge in range(0, cion_c):
        impdens = op_c.read_data_2d("DDLIMS", charge=charge+1, scaling=op_c.absfac)
        #ddlims_charges_c.append(impdens)
        #tot_dens += impdens
        zeff += impdens * (charge + 1)**2
    ddlims_charges_si = []
    for charge in range(0, cion_si):
        impdens = op_si.read_data_2d("DDLIMS", charge=charge+1, scaling=op_si.absfac)
        #ddlims_charges_si.append(impdens)
        #tot_dens += impdens
        zeff += impdens * (charge + 1)**2

    #for i in range(0, len(ne)):
    #    for charge in range(0, cion_c):
    #        zeff[i] += ddlims_charges_c[charge][i] * (charge + 1)
    #    for charge in range(0, cion_si):
    #        zeff[i] += ddlims_charges_si[charge][i] * (charge + 1)

    # Right now zeff is the total charge. Divide by tot_dens to get zeff.
    #zeff = zeff / tot_dens
    return zeff / ne

print("Calculating graphite Zeff...")
zeff_grap = zeff_grap(grap)
print("Calculating SiC Zeff...")
zeff_sic = zeff_sic(sic_c, sic_si)

# 2D plots.
grap.plot_contour_polygon("KTEBS", own_data=zeff_grap, cbar_label="Zeff (Graphite Walls)", vmin=1.0, vmax=1.5, cmap="inferno")
grap.plot_contour_polygon("KTEBS", own_data=zeff_sic, cbar_label="Zeff (SiC Walls)", vmin=1.0, vmax=1.5, cmap="inferno")
grap.plot_contour_polygon("KTEBS", own_data=zeff_grap-zeff_sic, cbar_label=r"$\mathdefault{Z_{eff}}$ Change", vmin=-0.1, vmax=0.1, cmap="coolwarm")
