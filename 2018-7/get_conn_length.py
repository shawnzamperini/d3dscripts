import scipy.interpolate as scinter
import numpy           as np
import get_ts          as ts
import time
import matplotlib.pyplot as plt
from scipy import spatial


def get_conn(shot, loc_time, loc_R, loc_Z, timeout=5, plot_it=True, grid_res=250, fig_ax1 = None):
    """
    Documentation.
    """
    # First get gfile (note this assumes remote connection).
    gfile = ts.load_gfile_mds(shot, loc_time)
    print("Gfile loaded.")

    # Get the relevant arrays into variables.

    # 2D array of wall coord. (R, Z).
    wall_RZ = gfile['wall']
    # Grid of R, Z values.
    Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
    # psin values at each (R, Z) location.
    psin_RZ = gfile['psiRZn']

    # R, Z of the magnetic axis.
    R_axis = gfile['RmAxis']
    Z_axis = gfile['ZmAxis']

    # Increase the resolution of the wall to a wall_res x wall_res grid.
    print('Increasing resolution between wall points...')
    wall_res = 100
    increased_wall_RZ = np.empty((0,2))
    for point in range(0, len(wall_RZ)-1):
        point1 = wall_RZ[point]
        point2 = wall_RZ[point+1]
        r1 = point1[0]
        z1 = point1[1]
        r2 = point2[0]
        z2 = point2[1]

        # If r1 = r2 then there's no distance between to interpolate to create the
        # new z's over. So instead get new r's from the z's.
        if r1 != r2:
            f_between         = scinter.interp1d([r1, r2], [z1, z2])
            new_wall_Rs       = np.linspace(r1, r2, wall_res)
            new_wall_Zs       = f_between(new_wall_Rs)
            more_wall_RZs     = np.column_stack((new_wall_Rs, new_wall_Zs))
            increased_wall_RZ = np.vstack((increased_wall_RZ, more_wall_RZs))
        else:
            f_between         = scinter.interp1d([z1, z2], [r1, r2])
            new_wall_Zs       = np.linspace(z1, z2, wall_res)
            new_wall_Rs       = f_between(new_wall_Zs)
            more_wall_RZs     = np.column_stack((new_wall_Rs, new_wall_Zs))
            increased_wall_RZ = np.vstack((increased_wall_RZ, more_wall_RZs))

    print("Creating interpolations and grids...")

    # Interpolation functions of psin(Z, R) and R(psin, Z).
    f_psin = scinter.Rbf(Rs, Zs, psin_RZ)
    f_R    = scinter.Rbf(psin_RZ, Zs, Rs, epsilon=0.00001)

    # Make a grid with a finer resolution.
    new_Rs = np.linspace(Rs.min(), Rs.max(), grid_res)
    new_Zs = np.linspace(Zs.min(), Zs.max(), grid_res)
    new_Rs, new_Zs = np.meshgrid(new_Rs, new_Zs)

    # Find the corresponding psins for this finer grid.
    new_psin_RZ      = f_psin(new_Rs, new_Zs)
    new_psin_RZ_flat = new_psin_RZ.flatten()

    # First find out what psin value loc_r and loc_z are on.
    loc_psin = f_psin(loc_R, loc_Z)

    # Create a 'KDTree' which when we give it a point, it returns the distance to
    # the closest point and the index of the point in the original array used
    # to create it, i.e. combined_RZ.
    combined_RZ = np.dstack([new_Rs.ravel(), new_Zs.ravel()])[0]
    mytree      = spatial.cKDTree(combined_RZ)

    # Get the grid point closest to the initial starting point.
    dist, idx_init = mytree.query([loc_R, loc_Z])
    print("Requested starting point: ({:.4f}, {:.4f})".format(loc_R, loc_Z))
    print("Closest grid point is:    ({:.4f}, {:.4f})".format(combined_RZ[idx_init][0], combined_RZ[idx_init][1]))

    # This index will be the same to get the corresponding psin in new_psin_RZ_flat.
    loc_psin = new_psin_RZ_flat[idx_init]
    print("Starting psin: {:.3f}".format(loc_psin))

    # Finding the connection length/distance to nearest limiter:
    #     1. There are 8 points surrounding every point. Find the one that has
    #          the value closest to loc_psin and move to it.
    #     2. Calculate the distance to this new point and add it to the
    #          connection length.
    #     3. If the new point is within wall_thresh of a wall point, break.
    #     4. Return to step 1 until you hit a wall.

    # Variables and arrays needed for the connection length loop. ITF.
    wall_thresh     = 0.01
    start_time      = time.time()
    current_R       = combined_RZ[idx_init][0]
    current_Z       = combined_RZ[idx_init][1]
    grid_spacing_R  = new_Rs[0][1] - new_Rs[0][0]
    grid_spacing_Z  = new_Zs[:,0][1] - new_Zs[:,0][0]
    conn_length_itf = 0
    closest_point   = None
    along_Rs_itf    = np.array([loc_R])
    along_Zs_itf    = np.array([loc_Z])
    start_up        = True
    dist = {}; idx = {}; tmp_psins = {}

    print("ITF: Finding distance to nearest limiter...")
    while True:
        # Get the 8 points surrounding the current point (l=lower/left, m=middle, u=upper, r=right).
        dist['ul'], idx['ul'] = mytree.query([current_R - grid_spacing_R, current_Z + grid_spacing_Z])
        dist['um'], idx['um'] = mytree.query([current_R                 , current_Z + grid_spacing_Z])
        dist['ur'], idx['ur'] = mytree.query([current_R + grid_spacing_R, current_Z + grid_spacing_Z])
        dist['ml'], idx['ml'] = mytree.query([current_R - grid_spacing_R, current_Z                 ])
        dist['mr'], idx['mr'] = mytree.query([current_R + grid_spacing_R, current_Z                 ])

        # Skip the bottom three on the first iteration since we know the ITF
        # will be the first point upwards, then set start_up to False later.
        if not start_up:
            dist['ll'], idx['ll'] = mytree.query([current_R - grid_spacing_R, current_Z - grid_spacing_Z])
            dist['lm'], idx['lm'] = mytree.query([current_R                 , current_Z - grid_spacing_Z])
            dist['lr'], idx['lr'] = mytree.query([current_R + grid_spacing_R, current_Z - grid_spacing_Z])

        # Get their corresponding psin value.
        tmp_psins['ul'] = new_psin_RZ_flat[idx['ul']]
        tmp_psins['um'] = new_psin_RZ_flat[idx['um']]
        tmp_psins['ur'] = new_psin_RZ_flat[idx['ur']]
        tmp_psins['ml'] = new_psin_RZ_flat[idx['ml']]
        tmp_psins['mr'] = new_psin_RZ_flat[idx['mr']]
        if not start_up:
            tmp_psins['ll'] = new_psin_RZ_flat[idx['ll']]
            tmp_psins['lm'] = new_psin_RZ_flat[idx['lm']]
            tmp_psins['lr'] = new_psin_RZ_flat[idx['lr']]

        # Don't want to move back to the previous point, so delete it for this
        # iteration. closest_point will still have last iteration's value.
        if closest_point is not None:
            if closest_point == 'ul':
                del tmp_psins['lr']
            elif closest_point == 'um':
                del tmp_psins['lm']
            elif closest_point == 'ur':
                del tmp_psins['ll']
            elif closest_point == 'ml':
                del tmp_psins['mr']
            elif closest_point == 'mr':
                del tmp_psins['ml']
            elif closest_point == 'll':
                del tmp_psins['ur']
            elif closest_point == 'lm':
                del tmp_psins['um']
            elif closest_point == 'lr':
                del tmp_psins['ul']
            else:
                print("How did you get here?")

        # Find out which point is the closest to loc_psin.
        tmp_psins     = {key: np.abs(tmp_psins[key] - loc_psin) for key in tmp_psins.keys()}
        closest_point = min(tmp_psins, key=tmp_psins.get)
        point         = combined_RZ[idx[closest_point]]
        how_far       = np.sqrt((current_R - point[0])**2 + (current_Z - point[1])**2)

        # Add to connection length.
        conn_length_itf  += how_far

        # Add points to bank for plotting later.
        along_Rs_itf = np.append(along_Rs_itf, point[0])
        along_Zs_itf = np.append(along_Zs_itf, point[1])

        # If within wall_thresh of wall, break.
        dist_from_wall = np.sqrt((point[0] - increased_wall_RZ[:,0])**2 + (point[1] - increased_wall_RZ[:,1])**2)
        if np.any(dist_from_wall < wall_thresh):
            print('Terminated on wall.')
            break

        # Timeout condition.
        if time.time() - start_time > timeout:
            print("Error: Timeout")
            break

        # Store for next iteration.
        current_R = point[0]
        current_Z = point[1]
        start_up = False

    # Variables and arrays needed for the connection length loop. OTF.
    wall_thresh     = 0.01
    start_time      = time.time()
    current_R       = combined_RZ[idx_init][0]
    current_Z       = combined_RZ[idx_init][1]
    grid_spacing_R  = new_Rs[0][1] - new_Rs[0][0]
    grid_spacing_Z  = new_Zs[:,0][1] - new_Zs[:,0][0]
    conn_length_otf = 0
    closest_point   = None
    along_Rs_otf    = np.array([loc_R])
    along_Zs_otf    = np.array([loc_Z])
    start_down      = True
    dist = {}; idx = {}; tmp_psins = {}

    print("OTF: Finding distance to nearest limiter...")
    while True:
        # Get the 8 points surrounding the current point (l=lower/left, m=middle, u=upper, r=right).
        if not start_down:
            dist['ul'], idx['ul'] = mytree.query([current_R - grid_spacing_R, current_Z + grid_spacing_Z])
            dist['um'], idx['um'] = mytree.query([current_R                 , current_Z + grid_spacing_Z])
            dist['ur'], idx['ur'] = mytree.query([current_R + grid_spacing_R, current_Z + grid_spacing_Z])
        dist['ml'], idx['ml'] = mytree.query([current_R - grid_spacing_R, current_Z                 ])
        dist['mr'], idx['mr'] = mytree.query([current_R + grid_spacing_R, current_Z                 ])
        dist['ll'], idx['ll'] = mytree.query([current_R - grid_spacing_R, current_Z - grid_spacing_Z])
        dist['lm'], idx['lm'] = mytree.query([current_R                 , current_Z - grid_spacing_Z])
        dist['lr'], idx['lr'] = mytree.query([current_R + grid_spacing_R, current_Z - grid_spacing_Z])

        # Get their corresponding psin value.
        if not start_down:
            tmp_psins['ul'] = new_psin_RZ_flat[idx['ul']]
            tmp_psins['um'] = new_psin_RZ_flat[idx['um']]
            tmp_psins['ur'] = new_psin_RZ_flat[idx['ur']]
        tmp_psins['ml'] = new_psin_RZ_flat[idx['ml']]
        tmp_psins['mr'] = new_psin_RZ_flat[idx['mr']]
        tmp_psins['ll'] = new_psin_RZ_flat[idx['ll']]
        tmp_psins['lm'] = new_psin_RZ_flat[idx['lm']]
        tmp_psins['lr'] = new_psin_RZ_flat[idx['lr']]

        # Don't want to move back to the previous point, so delete it for this
        # iteration. closest_point will still have last iterations value.
        if closest_point is not None:
            if closest_point == 'ul':
                del tmp_psins['lr']
            elif closest_point == 'um':
                del tmp_psins['lm']
            elif closest_point == 'ur':
                del tmp_psins['ll']
            elif closest_point == 'ml':
                del tmp_psins['mr']
            elif closest_point == 'mr':
                del tmp_psins['ml']
            elif closest_point == 'll':
                del tmp_psins['ur']
            elif closest_point == 'lm':
                del tmp_psins['um']
            elif closest_point == 'lr':
                del tmp_psins['ul']
            else:
                print("How did you get here?")

        # Find out which point is the closest to loc_psin.
        tmp_psins     = {key: np.abs(tmp_psins[key] - loc_psin) for key in tmp_psins.keys()}
        closest_point = min(tmp_psins, key=tmp_psins.get)
        point         = combined_RZ[idx[closest_point]]
        how_far       = np.sqrt((current_R - point[0])**2 + (current_Z - point[1])**2)

        # Add to connection length.
        conn_length_otf  += how_far
        # Add points to bank for plotting later.
        along_Rs_otf = np.append(along_Rs_otf, point[0])
        along_Zs_otf = np.append(along_Zs_otf, point[1])

        # If within wall_thresh of wall, break.
        dist_from_wall = np.sqrt((point[0] - increased_wall_RZ[:,0])**2 + (point[1] - increased_wall_RZ[:,1])**2)
        if np.any(dist_from_wall < wall_thresh):
            print('Terminated on wall.')
            break

        # Timeout condition.
        if time.time() - start_time > timeout:
            print("Error: Timeout")
            break

        current_R = point[0]
        current_Z = point[1]
        start_down = False

    if plot_it:
        if fig_ax1 is None:
            fig = plt.figure(figsize=(5,8))
            ax1 = fig.add_subplot(111)
        else:
            fig, ax1 = fig_ax1

        # Plot the wall first.
        ax1.plot(increased_wall_RZ[:,0], increased_wall_RZ[:,1], 'k')

        # Plot the starting point.
        ax1.plot(loc_R, loc_Z, 'k.', ms=10)

        # Plot the length along the field line.
        ax1.plot(along_Rs_itf, along_Zs_itf, 'r')

        # Plot the length along the field line.
        ax1.plot(along_Rs_otf, along_Zs_otf, 'r')

        fig.show()

    # Return the ITF and OTF connection length.
    return conn_length_itf, conn_length_otf, (fig, ax1)

def get_conn_list(shot, loc_time, loc_Rs, loc_Z, timeout=5, plot_it=True, grid_res=250):
    loc_Rs = np.array(loc_Rs)

    conns_itf = np.array([])
    conns_otf = np.array([])
    fig_ax1 = None

    for loc_R in loc_Rs:
        conn_itf, conn_otf, fig_ax1 = get_conn(shot, loc_time, loc_R, loc_Z, fig_ax1=fig_ax1)
        conns_itf = np.append(conns_itf, conn_itf)
        conns_otf = np.append(conns_otf, conn_otf)

    return conns_itf, conns_otf, fig_ax1
