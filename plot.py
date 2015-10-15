"""Short module to plot the output of GenPK"""

import numpy as np
import matplotlib.pyplot as plt
import re
import glob


def load_genpk(path,box, o_nu = 0):
    """Load a GenPk format power spectum, plotting the DM and the neutrinos (if present)
    Does not plot baryons."""
    #Load DM P(k)
    o_m = 0.3
    matpow=np.loadtxt(path)
    path_nu = re.sub("PK-DM-","PK-nu-",path)
    scale=2*math.pi/box
    #Adjust Fourier convention to match CAMB.
    simk=matpow[1:,0]*scale
    Pk=matpow[1:,1]/scale**3*(2*math.pi)**3
    return (simk,Pk)

def plot_genpk_power(matpow1, box,o_nu = 0, colour="blue"):
    """ Plot the matter power as output by gen-pk"""
    (k, Pk1)=load_genpk(matpow1,box, o_nu)
    #^2*2*!PI^2*2.4e-9*k*hub^3
    plt.ylabel("P(k) /(h-3 Mpc3)")
    plt.xlabel("k /(h Mpc-1)")
    plt.title("Power spectrum")
    plt.semilogx(k, Pk1, linestyle="-", color=colour)

