import sys
sys.path.append("/Users/zamperini/github/utk-fusion/tools")
import get_lp


lp_xl_path = "/Users/zamperini/Documents/d3d_work/184527/lp_184527.xlsx"
lp_xl_sheet = "Second Try"
lp_xl_xdata = "Psin.1"
lp_xl_ydata = "jsat (A/m2).1"
gauss_range= [0.998, 1.04]
lp_dict = get_lp.fit_conv_gauss(lp_xl_path=lp_xl_path, lp_xl_sheet=lp_xl_sheet,
  lp_xl_xdata=lp_xl_xdata, lp_xl_ydata=lp_xl_ydata, skiprows=1,
  gauss_range=gauss_range)

print("Psin")
for i in lp_dict["psin_fit"]:
    print(i)
print()
print("Y Data")
for i in lp_dict["y_fit"]:
    print(i)
