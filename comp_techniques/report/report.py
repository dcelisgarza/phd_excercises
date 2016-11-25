# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import math

# Expanded the sheath class from assignment 2. I couldn't find a way to do this with abstract classes.
class sheath:
    # Initialisation subroutine.
    def __init__(self, f0, vs = 1., mi = 1840., me = 1.):
        # Initialise f0 with the values of the starting conditions.
        # f0 = [phi_initial, E_initial]
        self.f0 = [f0[i] for i in range(len(f0))]    
        # vs = initial ion velocity
        self.vs  = vs
        # vs2 = vs^2
        self.vs2 = vs*vs
        # Exponential constant, mi and me take defaults if not given.
        self.cnt = np.sqrt(mi/(2.*np.pi*me))
        # Solutions.
        self.f = np.array([])
        # Currentn
        self.j = np.array([])
        
    # Checking initialised values.
    def show(self):
        print("Boundary conditions: ", self.f0)
        print("Initial ion velocity: ", self.vs)
        print("Initial ion velocity squared: ", self.vs2)
        print("Pre-Exponential constant: ", self.cnt)
    
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
    
    # Integrate.
    def integrate(self, x):
        self.f = sc.integrate.odeint(self, self.f0, x)
        
    # Calculate the current.
    def current(self):
        self.j = self.cnt * np.exp(self.f[:,0]) - 1.

class sheath_l(sheath):
    # Variable initialisation.
    def __init__(self, sheath, l):
        # Initialise f0 with the values of the starting conditions.
        # f0 = [phi_initial, E_initial, vi_initial]
        self.f0 = sheath.f0   
        # vs = initial ion velocity
        self.vs  = sheath.vs
        # vs2 = vs^2
        self.vs2 = sheath.vs2
        # Exponential constant, mi and me take defaults if not given.
        self.cnt = sheath.cnt
        # Variable length.
        self.l = l
        # Solutions.
        self.f = np.array([])
        # Currentn
        self.j = np.array([])
    
    def show_l(self):
        self.show()
        print("Ion collision length: ", self.l)
    
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

class sheath_r(sheath_l):
    # Variable initialisation.
    def __init__(self, sheath_l, r):
        # Initialise f0 with the values of the starting conditions.
        # f0 = [phi_initial, E_initial, vi_initial, ni_initial]
        self.f0 = sheath_l.f0
        # Exponential constant, mi and me take defaults if not given.
        self.cnt = sheath_l.cnt
        # Variable length.
        self.l = sheath_l.l
        # Recombination rate.
        self.r = r
        # Solutions.
        self.f = np.array([])
        # Currentn
        self.j = np.array([])
    
    def show_r(self):
        self.show_l()
        print("Recombination rate: ", self.r)
        
    # Subroutine called by the integrator.
    def __call__(self, f, x):
        # Electrostatic potential.
        phi = f[0]
        # Electric field.
        e   = f[1]
        # Ion velocity.
        vi  = f[2]
        # Ion number density.
        ni  = f[3]
        
        # Gradient of phi
        dphidx = -e
        # Ion velocity.
        dvidx  = e/vi - vi/self.l
        # dedx = -(d^2/dx^2) * phi = ni - ne
        dedx   = ni - np.exp(phi)
        # Change in ion number density.
        niovi = ni/vi
        dnidx  = self.r * niovi - niovi * dvidx
        return [dphidx, dedx, dvidx, dnidx]

def run():
    
    # Grid function.
    def grids(fig, alpha, minor = 0): 
        fig.grid(b = True, which = "major", linestyle = "--", alpha = alpha)
        if minor != 0:
            fig.grid(b = True, which = "minor", linestyle = "-.", alpha = alpha/4)
    
    # Close all figures.
    plt.close("all")
    # Matplotlib parameters to use TeX font and customise text size.
    plt.rc("text", usetex = True)
    plt.rc("font", family = "serif", size = 30)
    plt.rcParams['lines.linewidth'] = 4
    plt.rcParams["text.latex.preamble"] = r"\usepackage{bm}"
    #text.latex.preamble
    
    # Collision length parameter. (clp != 0, else we get 1/0)
    clp = eval(input("Collision length exponent (L = 10^x) = "))
    # Create figure.
    nplots = range(7)
    figarr = [plt.figure(i, figsize = (13,13/1.618)) for i in nplots]
    axarr  = [figarr[i].add_subplot(111) for i in nplots] 
    
    # Linspace for x.
    x = np.linspace(-40., 40., 10000)
    # Initialise empty array for vi @ wall.
    vi_wall = np.array([])
    # Logarithmic collision length array.
    lspace = np.array([clp,clp,1])
    l_arr = np.logspace(lspace[0], lspace[1], lspace[2])
    # Recombination rate:
    r_arr = np.array([0.])
    r_arr = np.append(r_arr, -np.logspace(math.copysign(lspace[0]+.15, -1), math.copysign(lspace[1]+.15, -1), lspace[2]))
    # Limits
    xlim = np.array([[0.,-0.], [0.,-0.], [0.,-0.]])
    xlimtmp = np.array([0.,-0.])
    xlimtmp2 = np.array([0.,-0.])
    ylim = np.array([[0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.]])
    ylimtmp = np.array([[0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.]])
    ylimtmp2 = np.array([[0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.]])
    ylimtmp3 = np.array([[0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.], [0.,-0.]])
    # Loop over l_arr
    for r in r_arr:
        for l in l_arr:
            # Initialise system with desired l.
            label1 = r"$L = 10^{%1i}, ~R = %1i$" %(np.log10(l), r)
            label2 = r"$L = 10^{%1i},~R = -10^{%1.2f}$" % (np.log10(l), math.copysign(np.log10(np.abs(r)),r))
            system = sheath([0., 0.001, 1., 1.])
            system = sheath_l(system, l)
            system = sheath_r(system, r)
            #system = sheath_l(sheath([0., 0.001, 1.], 1.), l)
            # Integrate
            system.integrate(x)
            # Calculate J.
            system.current()
            # Interpolate to find x when j = 0.
            x_wall = x - np.interp(0., system.j[::-1], x[::-1])
            ylimtmp[1:5,:] = [[np.interp(min(x_wall), x_wall, system.f[:,i]), np.interp(0., x_wall, system.f[:,i])] for i in range(4)]
            
            ylimtmp[6,:] = [np.interp(min(x_wall), x_wall, system.f[:,2]*system.f[:,3]), np.interp(0., x_wall, system.f[:,2]*system.f[:,3])]
            
            # Interpolate to find vi when x_wall = 0, vi_wall = vi(x_wall=0).
            # Append value to vi_wall array.
            vi_wall = np.append(vi_wall, np.interp(0., x_wall, system.f[:,2]))
            # Plot each curve.
            if r == 0:
                axarr[0].plot(x_wall, system.j, label = label1)
                [axarr[i+1].plot(x_wall, system.f[:,i], label = label1) for i in range(4)]
                axarr[5].plot(np.log(system.f[:,2]), np.log(system.f[:,3]), label = label1)
                axarr[6].plot(x_wall,system.f[:,2]*system.f[:,3], label = label1)
            else:
                axarr[0].plot(x_wall, system.j, label = label2)
                [axarr[i+1].plot(x_wall, system.f[:,i], label = label2) for i in range(4)]
                axarr[5].plot(np.log(system.f[:,2]), np.log(system.f[:,3]), label = label2)
                axarr[6].plot(x_wall,system.f[:,2]*system.f[:,3], label = label2)
        
        ylimtmp2[1:5,:] = [[min(ylim[i+1,:]), max(ylim[i+1,:])] for i in range(4)]
        ylimtmp3[1:5,:] = [[min(ylimtmp[i+1,:]), max(ylimtmp[i+1,:])] for i in range(4)]
        
        ylimtmp2[3,0] = ylim[3,0]
        ylimtmp3[3,0] = min(system.f[:,2])
        
        ylimtmp2[4,0] = ylim[4,0]
        ylimtmp3[4,0] = min(ylimtmp[4,:])
        ylimtmp2[4,1] = ylim[4,1]     
        ylimtmp3[4,1] = max(system.f[:,3])
        
        xlimtmp = [min(xlim[2,:]), max(xlim[2,:])]        
        xlimtmp2 = [min(np.log(system.f[:,2])), max(np.log(system.f[:,2]))]
        ylimtmp2[5,:] = [min(ylim[5,:]), max(ylim[5,:])]
        ylimtmp3[5,:] = [min(np.log(system.f[:,3])), max(np.log(system.f[:,3]))]
        
        ylimtmp2[6,:] = [min(ylim[6,:]), max(ylim[6,:])]
        ylimtmp3[6,:] = [min(ylimtmp[6,:]), max(ylimtmp[6,:])]
        
        ylimtmp2[6,0] = ylim[6,0]
        ylimtmp3[6,0] = min(ylimtmp[6,:])
        
        ylim[0,:] = np.array([0., max(system.j)])
        ylim[1:5,:] = [[min(ylimtmp2[i+1,0], ylimtmp3[i+1,0]), max(ylimtmp2[i+1,1], ylimtmp3[i+1,1])] for i in range(4)]
        
        if ylimtmp2[3,0] <= 0.:
            ylim[3,0] = max(ylimtmp2[3,0], ylimtmp3[3,0])
        else:
            ylim[3,0] = min(ylimtmp2[3,0], ylimtmp3[3,0])
        
        if ylimtmp2[4,0] > 0.:
            ylim[4,0] = min(ylimtmp2[4,0], ylimtmp3[4,0])
        else:
            ylim[4,0] = max(ylimtmp2[4,0], ylimtmp3[4,0])
        ylim[4,1] = max(ylimtmp2[4,1], ylimtmp3[4,1])
        
        ylim[6,:] = [min(ylimtmp2[6,0], ylimtmp3[6,0]), max(ylimtmp2[6,1], ylimtmp3[6,1])+0.01]
        if ylimtmp2[6,0] < 0.95:
            ylim[6,0] = max(ylimtmp2[6,0], ylimtmp3[6,0])
        else:
            ylim[6,0] = min(ylimtmp2[6,0], ylimtmp3[6,0])
        
        
        ylim[5,:] = [min(ylimtmp2[5,0], ylimtmp3[5,0]), max(ylimtmp2[5,1], ylimtmp3[5,1])]
        
        xlim[0,:] = np.array([min(min(xlim[0,:]), min(x_wall)) ,0.])
        xlim[1,:] = np.array([min(x_wall), 0.])
        xlim[2,:] = [min(xlimtmp[0], xlimtmp2[0]), max(xlimtmp[1], xlimtmp2[1])]
        
        # Add grids.  
        [grids(ax, 0.5) for ax in axarr]
        
        # Labels.
        [axarr[i].legend(loc=0) for i in nplots]
        [axarr[i].set_xlabel(r"$\bm{x}$, [Debye Lengths]") for i in range(5)]
        
        axarr[0].set_ylabel(r"$\bm{J}$, [Normalised]")
        axarr[1].set_ylabel(r"$\bm{\phi}$, [Normalised]")
        axarr[2].set_ylabel(r"$\bm{E}$, [Normalised]")
        axarr[3].set_ylabel(r"$\bm{v_{i}}$, [Normalised]")
        axarr[4].set_ylabel(r"$\bm{n_{i}}$, [Normalised]")
        
        axarr[5].set_xlabel(r"$\bm{\ln(v_{i})}$")
        axarr[5].set_ylabel(r"$\bm{\ln(n_{i})}$, [Normalised]")
        axarr[6].set_xlabel(r"$\bm{x}$, [Debye Lengths]")
        axarr[6].set_ylabel(r"$\bm{n_{i}v_{i}}$, [Normalised]")
        
        
        # Limits.
        axarr[0].set_ylim(ylim[0,:])
        [axarr[i+1].set_ylim(ylim[i+1,0], ylim[i+1,1]) for i in range(4)]
        [axarr[i].set_xlim(xlim[0,:]) for i in range(5)]
        axarr[5].set_xlim(xlim[2,:])
        axarr[5].set_ylim(ylim[5,:])
        axarr[6].set_xlim(xlim[1,:])
        axarr[6].set_ylim(ylim[6,:])
    
    #plt.tight_layout()
    # Save figure.
    [figarr[j].savefig(filename = "%1i_%il.eps" % (j, np.log10(l_arr[0])), format = "eps") for j in range(7)]
    #plt.savefig(filename = "collisions.eps", format = "eps")
    # Show plot.
    plt.show()

# Run code.
run()
