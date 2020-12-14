import matplotlib
import numpy as np
import matplotlib.pyplot as plt

#Initialization
xmin = -1 
xmax = 1    # x interval [xmin, xmax]
ymin = -0.5
ymax = 1.5  # y interval [ymin, ymax]
nx = 100    # number of points in [xmin, xmax]
ny = 100    # number of points in [ymin, ymax]
#Create x, y axis
x = np.linspace(xmin, xmax, nx)
y = np.linspace(ymin, ymax, ny)
#Create mesh
X, Y = np.meshgrid(x, y)
#Define fuction 
def f(x, y):
    return 11-2*x-20*y+21*x*x+10*y*y
#Evaluate functional values
Z = f(X, Y)
# Plot Contour
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
CS = plt.contour(X, Y, Z, 6, colors='k',)
plt.clabel(CS, fontsize=9, inline=1)
plt.title("Countour Plot of \n" 
          r"$m_k(p_1,p_2)=11-2p_1-20p_2+21p_1^2+10p_2^2$" 
          "\n of the quadratic model " 
          r"$ (4.2) $" " with " r"$x=(0,-1)$")
plt.xlabel(r"$p_1$")
plt.ylabel(r"$p_2$")
plt.show()