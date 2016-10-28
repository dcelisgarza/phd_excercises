# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 12:52:06 2016

@author: Daniel Celis Garza
"""

import matplotlib.pyplot   as plt
from scipy.integrate       import odeint
from numpy import pi       as np_pi
from numpy import exp      as np_exp
from numpy import sqrt     as np_sqrt
from numpy import linspace as np_linspace
from numpy import array    as np_array

# Close all figures.
plt.close("all")
# Matplotlib parameters to use TeX font and customise text size.
plt.rc("text", usetex = True)
plt.rc("font", family = "serif", size = 16)

# Define a general function for derivatives.
def derivs(f, x, params):
    # Electrostatic potential.
    phi = f[0]
    # Electric field.
    e   = f[1]
    # Calculate vi.
    vi  = np_sqrt(params[0]**2 - 2*phi)
    # Derivatives (d2phidx2 is actually de/dx which is the negative of d^2phi/dx^2).
    dphidx   = -e
    d2phidx2 = params[0]/vi - np_exp(phi)
    # Result of the function in the order f is given.
    return [dphidx, d2phidx2]

# x-array
x = np_linspace(0, 40, 100)
# Starting conditions for f_initial = [phi_initial, E_initial].
fi = np_array([0, 0.001])
# Vs = 1 (as an array so we can use the generalised derivs function in the future)
vs = np_array([1.0])
# Integrate.
y = odeint(derivs, fi, x, args = (vs,))

# Make two subplots for the assignment.
plt.subplots(2, 2, figsize=(10, 10))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

# Plot electrostatic potential and electric field.
ax1.plot(x,y[:, 0], label = r"Electrostatic Potential, $\phi$", linestyle = "-")
ax1.plot(x,y[:, 1], label = r"Electric Field, $E$", linestyle = "-")
ax1.set_xlabel(r"Debye Lengths")
ax1.set_ylabel(r"Normalised Potential \& Electric Field")
ax1.grid(b = True, which = "major", linestyle = "--", alpha = 0.6)
# Minor ticks option.
#ax1.grid(b = True, which = "minor", linestyle = "-.", alpha = 0.05)
#ax1.minorticks_on()
ax1.legend(loc = 3)

# Plot current.
mi = 1840
me = 1
j  = np_sqrt(mi/(2*np_pi*me)) * np_exp(y[:,0]) - 1

ax2.plot(x, j, label = r"Current, $J$")
ax2.set_xlabel(r"Debye Lengths")
ax2.set_ylabel(r"Normalised Current")
ax2.grid(b = True, which = "major", linestyle = "--", alpha = 0.6)
# Minor ticks option.
#ax2.grid(b = True, which = "minor", linestyle = "-.", alpha = 0.05)
#ax2.minorticks_on()
ax2.legend(loc = 0)

# Tighten layout.
plt.tight_layout()
# Show figure.
plt.show()
# Save figure.
plt.savefig(filename = "debye.eps", format = "eps")