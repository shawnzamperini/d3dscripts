# Helper script to read in a 3D wall file as a python dictionary, where each
# key is the toroidal angle from 0-359 degrees.

def read_wall(wall_path):

    # Read in the 3D wall coordinates.
    wall = {}
    with open(wall_path) as f:
        count = 0
        tor_angle = 0
        rs = []; zs = []
        for line in f:

            # Skip first three lines of header info and the first # of coordinates.
            if count < 4:
                count += 1
                continue

            # If we've hit a single number, this indicates the start of the next
            # set of coordinates, so save what we got and restart arrays.
            ls = line.split()
            if len(ls) == 1:

                # Add on the first point again so it completes the loop.
                rs.append(rs[0])
                zs.append(zs[0])

                wall[tor_angle] = [rs, zs]
                rs = []; zs = []
                tor_angle += 1

            # Save the pair of values in our coordinate lists.
            elif len(ls) == 2:
                rs.append(float(ls[0]))
                zs.append(float(ls[1]))
            else:
                print("Error! Line not recognized: {}".format(ls))
                break
            count += 1

    # Save the last entry since it won't have a number after it to trigger the
    # saving code.
    wall[tor_angle] = [rs, zs]

    # Add a 360 entry, which is just identical to 0 so things come full circle.
    wall[360] = wall[0]

    return wall
