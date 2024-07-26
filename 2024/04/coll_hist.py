import matplotlib.pyplot as plt
import pickle


path = "/Users/zamperini/flandir/coll_on/displ_dict_less_particles.pickle"
with open(path, "rb") as f:
    d = pickle.load(f)


