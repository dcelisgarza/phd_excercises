# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:47:14 2016

@author: dcg513
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

class sheath(object):
    
    # Initialisation subroutine.
    def __init__(self, f0, vs, l, mi = 1840., me = 1.):
        # Initialise f0 with the values of the starting conditions.
        # f0 = [phi_initial, E_initial, vi_initial]
        self.f0 = [f0[i] for i in range(len(f0))]
        # vs = initial ion velocity
        self.vs  = vs
        # vs2 = vs^2
        self.vs2 = vs*vs
        # l = mean collision length
        self.l = l
        # Exponential constant, mi and me take defaults if not given.
        self.cnt = np.sqrt(mi/(2.*np.pi*me))
    
    # Checking initialised values.
    def show(self):
        print(self.f0)
        print(self.vs)
        print(self.vs2)
        print(self.l)
        print(self.cnt)
    
    # Subroutine called by the integrator.
    def __call__(self, f, x):
        # Electrostatic potential.
        phi = f[0]
        # Electric field.
        e   = f[1]
        # Ion velocity.
        vi  = f[2]
        
        dphidx    = -e
        # dedx = -(d^2/dx^2) * phi
        dedx = self.vs/vi - np.exp(phi)
        dvidx = e/vi - vi/self.l
        return [dphidx, dedx, dvidx]
    
    # Calculate the current.
    def current(self, phi):
        return self.cnt * np.exp(phi) - 1

def run(system):
    
    # Grid function.
    def grids(alpha, minor = 0): 
        plt.grid(b = True, which = "major", linestyle = "--", alpha = alpha)
        if minor != 0:
            plt.grid(b = True, which = "minor", linestyle = "-.", alpha = alpha/4)
    
    # Close all figures.
    plt.close("all")
    # Matplotlib parameters to use TeX font and customise text size.
    plt.rc("text", usetex = True)
    plt.rc("font", family = "serif", size = 16)
    
    # Create figure.
    fig = plt.figure(1, figsize = (13,13/1.618))
    # Linspace for x.
    x = np.linspace(0., 100., 1000)
    # Initialise empty array for vi @ wall.
    vi_wall = np.array([])
    # Logarithmic collision length array.
    l_arr = np.logspace(-1,4,6)
    # Loop over l_arr
    for l in l_arr:
        # Initialise system with desired l.
        system = sheath([0., 0.001, 1.], 1., l)
        # Integrate
        f = sc.integrate.odeint(system, system.f0, x)
        # Calculate J.
        j = system.current(f[:,0])
        # Interpolate to find x when j = 0.
        x_wall = x - np.interp(0., j[::-1], x[::-1])
        # Interpolate to find vi when x_wall = 0.
        # Append value to vi_wall array.
        vi_wall = np.append(vi_wall, np.interp(0., x_wall, f[:,2]))
        # Plot each curve.
        plt.plot(x_wall, f[:,2], label = (r"$\hat{L} = %.1f$" % l))
        
    grids(0.5)
    plt.ylabel(r"Normalised Ion Velocity")
    plt.xlabel(r"Debye Lengths")
    plt.xlim([-20, 40]) ; plt.ylim([0, 5])
    plt.legend(loc=1)    
    
    ax = fig.add_axes([.125, .55, .15, .35])
    ax.plot(l_arr, vi_wall)
    ax.grid(alpha = 0.5)
    ax.set_ylabel(r"Ion Velocity @ Wall, $v_{i}(x = 0)$")
    ax.set_xlabel(r"Collision Length, L")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(filename = "collisions.eps", format = "eps")

system = sheath([0.,0., 0.], 0., 0.)
run(system)