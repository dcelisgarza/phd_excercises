import matplotlib.pyplot as plt
from numpy import linspace
from scipy.integrate import odeint

plt.close("all")

def lv(state, time, alpha, beta, gamma, delta):
    x = state[0] ; y = state[1]
    dxdt =  alpha * x - beta  * x * y
    dydt = -gamma * y + delta * x * y
    return [dxdt, dydt]


initial = [5,5]
t = linspace(0.,20,200)
result = odeint(lv,initial,t,args=(1.,1.,3.,1.))
plt.plot(t, result[:,0], label = "prey")
plt.plot(t, result[:,1], label = "predator")
plt.legend() ; plt.show()

#def getForce(x,y):
 #   N = len(x) ; Fx = zeros(N) ; Fy = zeros(N)
  #  for i in range(N):
   #     for j in range(N):
    #        if i == j: continue
     #       dx = x[i] - x[j] ; dy = y[i] - y[j]
        