# This is kind of just a scratch script to generate a PDF for input to DIVIMP
# option T29. It's pretty arbitrary since you need real data.
# It combines a skewed Gaussian outward radial velocoty PDF with an inward
# pinch Gaussian.
import numpy as np
from scipy.stats import skewnorm, norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Inputs
vmin         = 0
vmax         = 1000
vmean        = 400   # Mean of the Gaussian before applying the skew, see printout for true mean.
skew         = 1
blob_scale   = 1.0     # Multiplies the blob velocity PDF by this value. Hypothesis being that inertia of ions means not actually traveling at the blob velocity.
pinch        = -1
pinch_width  = 30
#pinch_weight = 1 - 1e-4   # When combining the PDFs how much total probability does the pinch get.
pinch_weight = 0.0

# Apply blob_scale.
vmin = vmin * blob_scale
vmax = vmax * blob_scale
vmean = vmean * blob_scale

# Skew Guassian outward radial velocity PDF.
vels = np.linspace(vmin, vmax, 50)
pdf1 = skewnorm.pdf(vels, skew, loc=vmean, scale=vels.max()/5)
print("Skewnorm mean: {:.2f}".format(skewnorm.mean(skew, loc=vmean, scale=vels.max()/5)))

# Normal Guassian pinch PDF.
pinches = np.linspace(pinch-pinch_width/2*2, pinch+pinch_width/2*2, 20)
pdf2 = norm.pdf(pinches, loc=pinch, scale=np.abs(pinch))

# Normalized for plotting individually.
pdf1_norm = (pdf1 - pdf1.min()) / (pdf1.max() - pdf1.min())
pdf2_norm = (pdf2 - pdf2.min()) / (pdf2.max() - pdf2.min())

# Create a common domain to sum the PDFs on.
vcom_min = min(vels.min(), pinches.min())
vcom_max = max(vels.max(), pinches.max())
vcom = np.linspace(vcom_min, vcom_max, 100)

# Create interpolatinos of each PDF, then interpolate onto common domain and
# sum accordingly.
f_pdf1 = interp1d(vels, pdf1, bounds_error=False, fill_value=0.0)
f_pdf2 = interp1d(pinches, pdf2, bounds_error=False, fill_value=0.0)
pdf1_com = f_pdf1(vcom)
pdf2_com = f_pdf2(vcom)
pdf_comb = pdf1_com * (1-pinch_weight) + pdf2_com * pinch_weight
pdf_comb_norm = (pdf_comb - pdf_comb.min()) / (pdf_comb.max() - pdf_comb.min())

# Print in input friendly format.
for i in range(0, len(vcom)):
    print("{:7.2f}  {:.2e}".format(vcom[i], pdf_comb_norm[i]))

# Plot them.
plt.rcParams['font.family'] = 'Century Gothic'
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))

ax1.plot(vels, pdf1_norm, label="Blob PDF", color="tab:red")
ax1.plot(pinches, pdf2_norm, label="Pinch PDF", color="tab:purple")
ax1.legend()
ax1.set_ylabel("Probability")
ax1.set_xlabel("Velocity (m/s)")
ax1.set_title("Separate PDFs")

ax2.plot(vcom, pdf_comb_norm, color="tab:red")
ax2.set_xlabel("Velocity (m/s)")
ax2.set_title("Combined PDF")

fig.tight_layout()
fig.show()
