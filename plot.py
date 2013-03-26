"""Short module to plot the output of GenPK"""

import numpy as np
import matplotlib.pyplot as plt
import re
import glob


def load_genpk(path,box, o_nu = 0):
    """Load a GenPk format power spectum."""
    #Load DM P(k)
    o_m = 0.3
    matpow=np.loadtxt(path)
    path_nu = re.sub("PK-DM-","PK-nu-",path)
    if glob.glob(path_nu):
        mp1a = np.loadtxt(path_nu)
        matpow_t = (mp1a*o_nu +matpow*(o_m - o_nu))/o_m
        ind = np.where(matpow_t/matpow > 2)
        matpow_t[ind] = matpow[ind]
        matpow = matpow_t
#         raise Exception
    scale=1.0/box
    #Adjust Fourier convention.
    simk=matpow[1:,0]*scale
    Pk=matpow[1:,1]/scale**3
    return (simk,Pk)

def plot_genpk_power(matpow1,matpow2, box,o_nu = 0, colour="blue"):
    """ Plot the matter power as output by gen-pk"""
    (k, Pk1)=load_genpk(matpow1,box, o_nu)
    (k,Pk2)=load_genpk(matpow2,box, o_nu)
    #^2*2*!PI^2*2.4e-9*k*hub^3
    plt.ylabel("P(k) /(h-3 Mpc3)")
    plt.xlabel("k /(h Mpc-1)")
    plt.title("Power spectrum change")
    plt.semilogx(k, Pk2/Pk1, linestyle="-", color=colour)


