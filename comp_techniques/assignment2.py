# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:42:14 2016

@author: dcg513
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

class sheath(object):
    
    def __init__(self, f0, vs, mi = 1840, me = 1):
        # Initialise f0 with the values of the starting conditions.
        self.f0 = [f0[i] for i in range(len(f0))]      
        self.vs  = vs
        self.vs2 = vs*vs
        self.cnt = np.sqrt(mi/(2*np.pi*me))
    
    def show(self):
        print(self.f0)
        print(self.vs)
        print(self.vs2)
        print(self.cnt)
    
    def __call__(self, f, x):
        phi = f[0]
        e   = f[1]
        vi  = np.sqrt(self.vs2 - 2*phi)
        dphidx    = -e
        nd2phidx2 = self.vs/vi - np.exp(phi)
        return [dphidx, nd2phidx2]
    
    def current(self, phi):
        return self.cnt * np.exp(phi) - 1

def run(system):
    
    def grids(alpha, minor = 1): # Grid function.
        plt.grid(b = True, which = "major", linestyle = "--", alpha = alpha)
        if minor == 1:
            plt.grid(b = True, which = "minor", linestyle = "-.", alpha = alpha/4)
    
    # Close all figures.
    plt.close("all")
    # Matplotlib parameters to use TeX font and customise text size.
    plt.rc("text", usetex = True)
    plt.rc("font", family = "serif", size = 16)
    
    plt.figure(1) # Create figure.
    x = np.linspace(0., 40., 100) # Create linspace.
    for vs in [1.0, 1.5, 2.0]: # Loop over vs
        system = sheath([0., 0.001], vs) # Initialise system with desired vs.
        f = sc.integrate.odeint(system, system.f0, x) # Integrate
        j  = system.cnt * np.exp(f[:,0]) - 1 # Calculate J
        x = x - np.interp(0., j[::-1], x[::-1]) # Interpolate to find x intersect.
        plt.plot(x, j, label = (r"$v_{s} = %.1f$" % system.vs)) # Plot each curve
    
    grids(0.5, 0.0)
    plt.ylabel(r"Normalised Current")
    plt.xlabel(r"Debye Lengths")
    plt.legend(loc=0)
    plt.tight_layout
    plt.tight_layout()
    plt.savefig(filename = "changevs.eps", format = "eps")


system = sheath([0.,0.001], 1.)
run(system)



