# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:42:14 2016

@author: dcg513
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

# Close all figures.
plt.close("all")
# Matplotlib parameters to use TeX font and customise text size.
plt.rc("text", usetex = True)
plt.rc("font", family = "serif", size = 16)

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
    
    def __call__(self, f, t):
        phi = f[0]
        e   = f[1]
        vi  = np.sqrt(self.vs2 - 2*phi)
        dphidx    = -e
        nd2phidx2 = self.vs/vi - np.exp(phi)
        return [dphidx, nd2phidx2]
    
    def current(self, phi):
        return self.cnt * np.exp(phi) - 1

def run(system):
    
    def 
    def ax(ax, position):
        ax = plt.subplots(position)
        return ax
    
    def label(system, i):
        if i == 0:            
            return r"Electrostatic Potential, $\phi$"
        elif i == 1:
            return r"Electric Field, $E$"
        elif i == 3:
            return r"V_{s} = %.1f" % (system.vs)
        elif i == 4:
            return r"Current, $J$"
        else:
            return r"V_{s} = %.1f" % (system.vs)
    
    def axis_label(ax, i):
        if i == 0:
            return ax.set_xlabel(r"Debye Lengths")
        elif i == 1:
            return ax.set_ylabel(r"Normalised Potential \& Electric Field")
        else:
            return ax.set_ylabel(r"Normalised Current")
    
    def grids(ax, minor = True, alpha):
        ax.grid(b = True, which = "major", linestyle = "--", alpha = alpha)
        if minor == True:
            ax.grid(b = True, which = "minor", linestyle = "-.", alpha = alpha/2)
        

    x = np.linspace(0., 40., 100)
    f = sc.integrate.odeint(system, system.f0, x)
    plt.subplots(2, 2, figsize=(10, 10))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    
    for i in range(len(f[1])):
        ax1.plot(x, f[:,i], label = label(system, i))
    for i in [0, 1]:
        axis_label(i)
    
    ax1.grid(b = True, which = "major", linestyle = "--", alpha = 0.6)
    ax1.grid(b = True, which = "minor", linestyle = "-.", alpha = 0.3)
    ax1.minorticks_on()
    ax1.legend(loc=0)

test = sheath([0.,0.001], 1.)
test.show()
run(test)

    
# Define a general function for derivatives.
def derivs(f, x, params):
    # Electrostatic potential.
    phi = f[0]
    # Electric field.
    e   = f[1]
    # Calculate vi.
    vi  = np.sqrt(params[0]**2 - 2*phi)
    # Derivatives (d2phidx2 is actually de/dx which is the negative of d^2phi/dx^2).
    dphidx   = -e
    d2phidx2 = params[0]/vi - np.exp(phi)
    # Result of the function in the order f is given.
    return [dphidx, d2phidx2]

# x-array
x = np.linspace(0, 40, 100)
# Starting conditions for f_initial = [phi_initial, E_initial].
fi = np.array([0, 0.001])
# Vs = 1 (as an array so we can use the generalised derivs function in the future)
vs = np.array([1.0])
# Integrate.
y = sc.integrate.odeint(derivs, fi, x, args = (vs,))

# Make two subplots for the assignment.
plt.subplots(2, 2, figsize=(10, 10))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

# Plot electrostatic potential and electric field.
ax1.plot(x,y[:, 0], label = r"Electrostatic Potential, $\phi$", linestyle = "-")
ax1.plot(x,y[:, 1], label = r"Electric Field, $E$", linestyle = "--")
ax1.set_xlabel(r"Debye Lengths")
ax1.set_ylabel(r"Normalised Potential \& Electric Field")
ax1.grid(b = True, which = "major", linestyle = "--", alpha = 0.6)
# Minor ticks option.
ax1.grid(b = True, which = "minor", linestyle = "-.", alpha = 0.3)
ax1.minorticks_on()
ax1.legend(loc = 3)

# Plot current.current
mi = 1840
me = 1
j  = np.sqrt(mi/(2*np.pi*me)) * np.exp(y[:,0]) - 1

ax2.plot(x, j, label = r"Current, $J$")
ax2.set_xlabel(r"Debye Lengths")
ax2.set_ylabel(r"Normalised Current")
ax2.grid(b = True, which = "major", linestyle = "--", alpha = 0.6)
# Minor ticks option.
ax2.grid(b = True, which = "minor", linestyle = "-.", alpha = 0.3)
ax2.minorticks_on()
ax2.legend(loc = 0)

# Tighten layout.
plt.tight_layout()
# Show figure.

# Save figure.
plt.savefig(filename = "debye.eps", format = "eps")