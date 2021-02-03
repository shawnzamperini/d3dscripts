import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


xl_path = "/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/Slides, Sheets and Documents/2020/11/dissertation_data.xlsx"

fontsize = 16
ms = 18

def lambdas():

    df = pd.read_excel(xl_path, sheet_name="lambdas")

    x_rev = df["# of lambdas"][:5]
    y_rev = df["ITF/OTF Total"][:5]
    yerr_rev = df["ITF/OTF Error"][:5]
    x_for = df["# of lambdas"][5:]
    y_for = df["ITF/OTF Total"][5:]
    yerr_for = df["ITF/OTF Error"][5:]

    fig, ax = plt.subplots()
    ax.errorbar(x=x_rev, y=y_rev, color="tab:purple", label=r"Bx$\nabla$B$\uparrow$",   fmt='.', ms=ms, mec='k', mew=1, capsize=2, ecolor='k')
    ax.errorbar(x=x_for, y=y_for, color="tab:red",    label=r"Bx$\nabla$B$\downarrow$", fmt='.', ms=ms, mec='k', mew=1, capsize=2, ecolor='k')
    ax.set_xlabel("# of " + r'$\mathrm{\lambda_{ne}}$' + "'s from separatrix", fontsize=fontsize)
    ax.set_ylabel("ITF/OTF Total", fontsize=fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=fontsize)
    #ax.set_yscale("log")
    fig.tight_layout()
    fig.show()

    return df
