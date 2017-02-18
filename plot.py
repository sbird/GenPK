"""Short module to plot the output of GenPK"""

import numpy as np
import matplotlib.pyplot as plt
import re
import math


def load_genpk(path,box):
    """Load a GenPk format power spectum, plotting the DM and the neutrinos (if present)
    Does not plot baryons."""
    #Load DM P(k)
    matpow=np.loadtxt(path)
    scale=2*math.pi/box
    #Adjust Fourier convention to match CAMB.
    simk=matpow[1:,0]*scale
    Pk=matpow[1:,1]/scale**3*(2*math.pi)**3
    return (simk,Pk)

def plot_genpk_power(matpow1, box,color=None, ls="-"):
    """ Plot the matter power as output by gen-pk"""
    (k, Pk1)=load_genpk(matpow1,box)
    #^2*2*!PI^2*2.4e-9*k*hub^3
    plt.ylabel("P(k) /(h-3 Mpc3)")
    plt.xlabel("k /(h Mpc-1)")
    plt.title("Power spectrum")
    plt.loglog(k, Pk1, linestyle=ls, color=color)

def get_camb_power(matpow):
    """Plot the power spectrum from CAMB
    (or anything else where no changes are needed)"""
    data = np.loadtxt(matpow)
    return (data[:,0], data[:,1])
