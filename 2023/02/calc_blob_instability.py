# Using BlobbyFarSOL to calculate the neutral densities, calculate the final unknown in Theiler's thesis Eq. 5.1.13,
# the secondary linear growth rate that can limit its motion.
import sys
sys.path.append("../../2022/12")
import BlobbyFarSOL
import matplotlib.pyplot as plt


# Inputs
shot = 190484

bfs = BlobbyFarSOL.BlobbyFarSOL()