# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:42:14 2016

@author: dcg513
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# Sheath object for easy generalisation.
class sheath(object):
    
    # Initialisation subroutine.
    def __init__(self, f0, vs, mi = 1840., me = 1.):
        # Initialise f0 with the values of the starting conditions.
        # f0 = [phi_initial, E_initial]
        self.f0 = [f0[i] for i in range(len(f0))]    
        # vs = initial ion velocity
        self.vs  = vs
        # vs2 = vs^2
        self.vs2 = vs*vs
        # Exponential constant, mi and me take defaults if not given.
        self.cnt = np.sqrt(mi/(2.*np.pi*me))
    
    # Checking initialised values.
    def show(self):
        print(self.f0)
        print(self.vs)
        print(self.vs2)
        print(self.cnt)
    
    # Subroutine called by the integrator.
    def __call__(self, f, x):
        # Electrostatic potential.
        phi = f[0]
        # Electric field.
        e   = f[1]
        # Ion velocity.
        vi  = np.sqrt(self.vs2 - 2.*phi)
        
        dphidx = -e
        # dedx = -(d^2/dx^2) * phi
        dedx   = self.vs/vi - np.exp(phi)
        return [dphidx, dedx]
    
    # Calculate the current.
    def current(self, phi):
        return self.cnt * np.exp(phi) - 1.

def run(system):
    
    # Grid function.
    def grids(alpha, minor = 0): 
        plt.grid(b = True, which = "major", linestyle = "--", alpha = alpha)
        if minor != 0:
            plt.grid(b = True, which = "minor", linestyle = "-.", alpha = alpha/4.)
    
    # Close all figures.
    plt.close("all")
    # Matplotlib parameters to use TeX font and customise text size.
    plt.rc("text", usetex = True)
    plt.rc("font", family = "serif", size = 16)
    
    # Create figure.
    plt.figure(1) 
    # Linspace for x.
    x = np.linspace(0., 40., 100) 
    # vs array.
    vs_arr = np.linspace(1.,2.,3) 
    # Loop over vs_arr.
    for vs in vs_arr:
        # Initialise system with desired vs.
        system = sheath([0., 0.001], vs) 
        # Integrate.
        f = sc.integrate.odeint(system, system.f0, x)
        # Calculate J.
        j  = system.current(f[:,0])
        # Interpolate to find x when J = 0.
        x_wall = x - np.interp(0., j[::-1], x[::-1])
        # Plot each curve.
        plt.plot(x_wall, j, label = (r"$\hat{v_{s}} = %.1f$" % system.vs)) 
    
    # Add grids.
    grids(0.5)
    # Labels.
    plt.ylabel(r"$J$, [Normalised]")
    plt.xlabel(r"$x_{wall}$, [Debye Lengths]")
    # Legend position.
    plt.legend(loc=0)
    # Tight layout.
    plt.tight_layout()
    # Save figure.
    plt.savefig(filename = "changevs.eps", format = "eps")
    # Show figure.
    plt.show()

# Initialise system with any values, they're changed inside run.
system = sheath([0.,0.], 0.)
# Run code.
run(system)



