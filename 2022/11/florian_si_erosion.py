import pandas as pd
import sys
sys.path.append("/Users/zamperini/github/d3dscripts/2022/08/")
import EngelhardtModel


# Measured values from Florian's paper.
dflux = 2.1e18
siflux5 = 1.75e15
siflux15 = 2.98e15
eimpact5 = 65
eimpact15 = 70
fc = 0.02
include_carbon = True

if not include_carbon:
    fc = 0.0

# Load yields. We just want the Silicon ones really.
em = EngelhardtModel.EngelhardtModel()
em.load_mm_model()

total_yield5 = em.Y_D_Si(eimpact5) + fc*em.Y_C_Si(eimpact5) + em.Y_D_Si_ch300(eimpact5)
total_yield15 = em.Y_D_Si(eimpact15) + fc*em.Y_C_Si(eimpact15) + em.Y_D_Si_ch300(eimpact15)

gross5 = dflux * total_yield5
gross15 = dflux * total_yield15

print("Prompt Redeposition Fraction")
print(" 5Hz:  {:.3f}".format(1-(siflux5/gross5)))
print(" 15Hz: {:.3f}".format(1-(siflux15/gross15)))
