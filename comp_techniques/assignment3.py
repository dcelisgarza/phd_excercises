# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:47:14 2016

@author: dcg513
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

class sheath(object):
    
    def __init__(self, f0, vs, l, mi = 1840, me = 1):
        # Initialise f0 with the values of the starting conditions.
        self.f0 = [f0[i] for i in range(len(f0))]      
        self.vs  = vs
        self.vs2 = vs*vs
        self.l = l
        self.cnt = np.sqrt(mi/(2*np.pi*me))
    
    def show(self):
        print(self.f0)
        print(self.vs)
        print(self.vs2)
        print(self.l)
        print(self.cnt)
    
    def __call__(self, f, x):
        phi = f[0]
        e   = f[1]
        vi  = f[2]
        dphidx    = -e
        nd2phidx2 = self.vs/vi - np.exp(phi)
        dvidx = e/vi - vi/self.l
        return [dphidx, nd2phidx2, dvidx]
    
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
    
    viwall = np.array([])
    L = np.logspace(-1,4,6)
    fig = plt.figure(1, figsize = (13,13/1.618)) # Create figure.
    x = np.linspace(0., 100., 1000) # Create linspace.
    for l in L: # Loop over vs
        system = sheath([0., 0.001, 1.], 1., l) # Initialise system with desired vs.
        f = sc.integrate.odeint(system, system.f0, x) # Integrate
        j = system.current(f[:,0])
        xwall = x - np.interp(0., j[::-1], x[::-1]) # Interpolate to find x when j = 0.
        viwall = np.append(viwall, np.interp(0., xwall, f[:,2])) # Interpolate to find vi when xwall = 0.
        plt.plot(xwall, f[:,2], label = (r"$\hat{L} = %.1f$" % l)) # Plot each curve
        
    grids(0.5, 0.0)
    plt.ylabel(r"Normalised Ion Velocity")
    plt.xlabel(r"Debye Lengths")
    plt.xlim([-20, 40]) ; plt.ylim([0, 5])
    plt.legend(loc=1)    
    
    ax = fig.add_axes([.125, .55, .15, .35])
    ax.plot(L, viwall)
    ax.grid(alpha = 0.5)
    ax.set_ylabel(r"Ion Velocity @ Wall, $v_{i}(x = 0)$")
    ax.set_xlabel(r"Collision Length, L")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(filename = "collisions.eps", format = "eps")

system = sheath([0.,0.001, 1.], 1., 0.1)
run(system)