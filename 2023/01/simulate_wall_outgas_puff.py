# This script takes an OEDGE grid and divides a total amount of wall outgassing among a number of point sources
# along each wall element. It's a bit of a hack job, but there isn't much else we can do about it.
import oedge_plots
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

# Load run. Doesn't matter which since we just want the wall coordinates.
# ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/190484/d3d-190484-bkg-002.nc"
ncpath = "/Users/zamperini/Documents/d3d_work/divimp_files/167196/d3d-167196-bg-shifted-ring-entry-12.nc"
op = oedge_plots.OedgePlots(ncpath)
r0 = float(op.nc["R0"][:].data)
z0 = float(op.nc["Z0"][:].data)

# Pull out the vessel coordinates. This contains the R, Z coordinates of all the vertices.
keep_idx = np.where(np.logical_and(op.rvesm[0] != 0, op.zvesm[0] != 0))
rvesm_original = np.append(op.rvesm[0][keep_idx], op.rvesm[0][keep_idx][0])
zvesm_original = np.append(op.zvesm[0][keep_idx], op.zvesm[0][keep_idx][0])
vessel = Polygon(zip(rvesm_original, zvesm_original))


# Stolen from https://stackoverflow.com/questions/49558464/shrink-polygon-using-corner-coordinates
def shrink_or_swell_shapely_polygon(my_polygon, factor=0.10, swell=False):
    """ returns the shapely polygon which is smaller or bigger by passed factor.
        If swell = True , then it returns bigger polygon, else smaller """

    from shapely import geometry
    xs = list(my_polygon.exterior.coords.xy[0])
    ys = list(my_polygon.exterior.coords.xy[1])
    x_center = 0.5 * min(xs) + 0.5 * max(xs)
    y_center = 0.5 * min(ys) + 0.5 * max(ys)
    min_corner = geometry.Point(min(xs), min(ys))
    max_corner = geometry.Point(max(xs), max(ys))
    center = geometry.Point(x_center, y_center)
    shrink_distance = center.distance(min_corner) * factor

    if swell:
        my_polygon_resized = my_polygon.buffer(shrink_distance)  # expand
    else:
        my_polygon_resized = my_polygon.buffer(-shrink_distance)  # shrink

    # Visualize for debugging
    x, y = my_polygon.exterior.xy
    plt.plot(x, y)
    x, y = my_polygon_resized.exterior.xy
    plt.plot(x, y)
    # To not let the image be distorted along the axis
    plt.axis('equal')
    plt.show()

    return my_polygon_resized


# What I'm doing here is an absolute hack job that I'm not taking the time to fully understand but it does the job. I
# pass in the Polygon and then I shrink it just a small amount so that the points are slightly removed from the wall.
# The idea is to help get rid of some errors when inputting the data for EIRENE.
vessel_shrunk = shrink_or_swell_shapely_polygon(vessel, factor=0.01)
rvesm = np.array(list(vessel_shrunk.exterior.coords))[:, 0]
zvesm = np.array(list(vessel_shrunk.exterior.coords))[:, 1]

# The resolution of the points, or how far space along the wall each point source is. Due to weird bug, EIRENE caps
# at 34 points here, so make sure you're under that?
res = 0.23

# Source strength in Torr-L/s. We then convert to Amps as needed by EIRENE. WALL_RATE for 190484 is around 20 Torr-L/s.
f_tlps = 50
f_sccm = f_tlps * 78.9
f_pps = f_sccm * 4.477962 * 10 ** 17
f_amps = f_pps * 1.6022e-19

s = 0.0
step = res
source_data = {"r": [], "z": [], "vecx": [], "vecy": [], "numpoints": 0}
for i in range(0, len(rvesm)):

    # Calculate length of segment we're on.
    if i == len(rvesm) - 1:
        seg_len = np.sqrt(np.square(rvesm[0] - rvesm[i]) + np.square(zvesm[0] - zvesm[i]))

        # Catch vertical lines, assigning correct angle depending on if  we're going up or down.
        if (rvesm[0] - rvesm[i]) == 0:
            if zvesm[0] > zvesm[i]:
                seg_ang = np.pi / 2
            else:
                seg_ang = -np.pi / 2
        else:
            seg_ang = np.arctan2((zvesm[0] - zvesm[i]), (rvesm[0] - rvesm[i]))
    else:
        seg_len = np.sqrt(np.square(rvesm[i + 1] - rvesm[i]) + np.square(zvesm[i + 1] - zvesm[i]))
        if (rvesm[i + 1] - rvesm[i]) == 0:
            if zvesm[i + 1] > zvesm[i]:
                seg_ang = np.pi / 2
            else:
                seg_ang = -np.pi / 2
        else:
            seg_ang = np.arctan2((zvesm[i + 1] - zvesm[i]), (rvesm[i + 1] - rvesm[i]))

    # If the segment is too short, just subtract that segment length from our current step size
    # and move onto the next segment.
    if step > seg_len:
        step -= seg_len
        continue

    seg_s = 0.0
    r = rvesm[i]
    z = zvesm[i]
    while True:

        # If we ran out of segment, make it so the next step size is shorter by the amount of leftover distance.
        if (seg_len - seg_s) < step:
            step -= (seg_len - seg_s)
            break

        # Take a step along the surface and calculate the R, Z where you end up.
        r += step * np.cos(seg_ang)
        z += step * np.sin(seg_ang)

        # Add a point source at this location.
        source_data["r"].append(r)
        source_data["z"].append(z)

        # Calculate the vector value for this point. We choose it so the source is directed towards the plasma center.
        # If we use the comment in the input file that vecx=1.0 is inwards, and vecy=1.0 is downwards, we're good. This
        # is a strange convention I don't fully understand though.
        vecx = r - r0
        vecy = z - z0

        # Normalize the vectors so that one of them is 1.0 (and the other less than 1.0). Unnecessary since EIRENE
        # normalizes them, but makes input cleaner.
        maxval = max(np.abs(vecx), np.abs(vecy))
        source_data["vecx"].append(vecx / maxval)
        source_data["vecy"].append(vecy / maxval)

        # Update variables for next iteration.
        seg_s += step
        s += step
        step = res
        source_data["numpoints"] += 1

print("Created {} point sources.".format(source_data["numpoints"]))

# Put into a textfile the output in a format that is friendly for copy/pasting into the DIVIMP input file.
puff_str = f_amps / source_data["numpoints"]
puff_ev = 0.026
with open("wall_outgass_puff.dat", "w") as f:
    for i in range(0, source_data["numpoints"]):
        f.write("3.0  -90000  {:.2e}  1.0  2  1  {:.3f}  1.0  90.0  {:.3f}  {:.3f}  0.000  {:.2f}  {:.2f}  0.0  '{}'\n"
                .format(puff_str, puff_ev, source_data["r"][i], source_data["z"][i], source_data["vecx"][i],
                        source_data["vecy"][i], "wall outgassing {}".format(i)))

# A plot of the output.
fig, ax = plt.subplots()
ax.plot(rvesm_original, zvesm_original, color="k")

ax.scatter(source_data["r"], source_data["z"], s=15, color="r")

ax.set_aspect("equal")
ax.set_xlabel("R (m)")
ax.set_ylabel("Z (m)")
fig.tight_layout()
fig.show()
