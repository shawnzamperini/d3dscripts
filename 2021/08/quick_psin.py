import numpy  as np
from scipy.interpolate import interp1d, interp2d
import pickle

z = np.arange(-0.85, -1.255, -0.005)
r = np.full(len(z), 1.4859)

gfile_path = "/Users/zamperini/Google Drive/My Drive/Research/Data/rcp_data/gfile_data_for_rcp/" + \
  "{}".format("187107_1620")

with open(gfile_path, "rb") as f:
    gfile = pickle.load(f)

f_psin = interp2d(gfile["R"], gfile["Z"], gfile["Psin"])

rcp_psin = f_psin(r, z)[:, 0]
