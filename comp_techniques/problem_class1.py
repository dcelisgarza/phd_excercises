# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:19:01 2016

@author: dcg513
"""

#Import useful routines from modules
from scipy.integrate import odeint
from numpy import linspace, exp, sin, cos, array
import matplotlib.pyplot as plt

plt.close("all")
plt.rc("text", usetex = True)
plt.rc("font", family = "serif", size = 16)

#Here we define our function which returns df/dt = -a*f
#Note we can pass in a, but it defaults to 10
def dfdt(curF, curT, params):
    #We don’t do anything with curT
    return -params[0]*sin(curF) - params[1]*cos(curT)
#Now define the times at which we want to know the result
time=linspace(0,10,400)
#Set the initial condition
f0=10.

params = array([30.,50.])
result = odeint(dfdt,f0,time,args=(params,))
plt.figure(1)
plt.plot(time, result, label = r"$\frac{\mathrm{d}f}{\mathrm{d}t} = -" + str(params[0]) + " \sin(f) -" +str(params[1]) + " \cos(t)$")
plt.xlabel("Time"); plt.ylabel(r"$f$")
plt.legend(loc=0) ; plt.show()

def spring(f, x, params):
    # params = [k, m]
    # Assign x.
    x = f[0]
    # Assign v.
    v = f[1]
    # Construct derivatives.
    dxdt   = v
    d2xdt2 = -params[0]/params[1] * x
    # Return derivatives.
    return [dxdt, d2xdt2]

params = array([0.5,0.01])
f = array([0.4,0.])
time = linspace(0,10,1000)
result = odeint(spring,f,time,args=(params,))
plt.figure(2)
plt.plot(time, result[:,0], label = r"$x(t)$.", linestyle = "-")
plt.plot(time, result[:,1], label = r"$v(t)$.", linestyle = "--")
plt.legend(loc=0)
plt.xlabel(r"Time, s"); plt.ylabel(r"Value")
plt.figure(3)
plt.plot(time, params[1]*result[:,1]**2/2, label = r"$K(t)$", linestyle = "--")
plt.plot(time, params[0]*result[:,0]**2/2, label = r"$V(t)$", linestyle = "-.")
plt.plot(time, params[0]*result[:,0]**2/2 + params[1]*result[:,1]**2/2, label = r"$E(t)$", linestyle = "-")
plt.legend(loc=0)
plt.xlabel(r"Time, s"); plt.ylabel(r"Energy")

def dampspring(f, x, params):
    # params = [k, m]
    # Assign x.
    x = f[0]
    # Assign v.
    v = f[1]
    # Construct derivatives.
    dxdt   = v
    d2xdt2 = -params[0]/params[1] * x - params[2]/params[1]*v
    # Return derivatives.
    return [dxdt, d2xdt2]

params = array([0.5,0.01,0.005])
f = array([0.4,0.])
time = linspace(0,10,1000)
result = odeint(dampspring,f,time,args=(params,))
plt.figure(4)
plt.plot(time, result[:,0], label = r"$x(t)$.", linestyle = "-")
plt.plot(time, result[:,1], label = r"$v(t)$.", linestyle = "--")
plt.legend(loc=0)
plt.xlabel(r"Time, s"); plt.ylabel(r"Value")
plt.figure(5)
plt.plot(time, params[1]*result[:,1]**2/2, label = r"$K(t)$", linestyle = "--")
plt.plot(time, params[0]*result[:,0]**2/2, label = r"$V(t)$", linestyle = "-.")
plt.plot(time, params[0]*result[:,0]**2/2 + params[1]*result[:,1]**2/2, label = r"$E(t)$", linestyle = "-")
plt.legend(loc=0)
plt.xlabel(r"Time, s"); plt.ylabel(r"Energy")

    