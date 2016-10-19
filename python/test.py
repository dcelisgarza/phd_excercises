import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.stats import norm
plt.close("all")

# Plot between -10 and 10 with .001 steps.
x_axis = np.arange(-10, 10, 0.001)
# Mean = 0, SD = 2.
plt.plot(x_axis, norm.pdf(x_axis,0,2))

plt.show()


a = np.random.rand()

b = np.random.uniform()
c = np.random.exponential()

print(a,b,c)
