import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

###### Initialization #########
# This is the code for Part I: x = (0,-1), 
# i.e. x1 = 0, x2 = -1 Which is the starting point
x0 = 0
y0 = -1
gs = 3  # Grid Size
ng = 500    # Grid steps
Dmax = 2    # Maximum TR radius
nD = 200 # Number of steps for Trust radius
nth = 200 # Angle steps
# Plot the initial point

#Create x, y axis
xmin = x0-gs
xmax = x0+gs
ymin = y0-gs
ymax = y0+gs
x = np.linspace(xmin, xmax, ng)
y = np.linspace(ymin, ymax, ng)
#Create mesh
X, Y = np.meshgrid(x, y)
#Define fuction 
def f(x, y):
    return 10*np.multiply(((y-np.multiply(x,x))),((y-np.multiply(x,x)))) + np.multiply(1-x,1-x)
#Evaluate functional values
Z = f(X, Y)
f0 = 10*(y0-x0**2)**2 + (1-x0)**2;   # initial function values
g0 = np.matrix([[-40*x0*(y0-x0**2)-2*(1-x0)], [20*(y0-x0**2)]])
B0 = np.matrix([[120*x0**2-40*y0+2, -40*x0], [-40*x0, 20]])
M = f0 + g0.item(0)*(X-x0)+g0.item(1)*(Y-y0) + 0.5*B0.item((0,0))*np.multiply((X-x0),(X-x0)) + B0.item((0,1))*np.multiply((X-x0),(Y-y0)) + 0.5*B0.item((1,1))*np.multiply((Y-y0),(Y-y0))
# Plot Contour
plt.plot([x0], [y0], 'g*', markersize=15)
x0y0_legend = mlines.Line2D([], [], color='w', marker='*', markerfacecolor='g',
                          markersize=15, label='$(x_0,y_0)$')

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
CS = plt.contour(X, Y, Z, 30, colors='b', )
plt.clabel(CS, fontsize=9, inline=1)
f_legend = mlines.Line2D([], [], color='b', label='$f(x,y)$')

CS = plt.contour(X, Y, Z,levels = [f0],
                 colors=('k',),linestyles=('-',),linewidths=(2,))
plt.clabel(CS, fontsize=9, inline=1)
f_legend0 = mlines.Line2D([], [], color='k', label='$f(x_0,y_0)$')

CS = plt.contour(X, Y, M, 8, colors='r',)
plt.clabel(CS, fontsize=9, inline=1)
M_legend = mlines.Line2D([], [], color='r', label='$m_k(x,y)$')

CS = plt.contour(X, Y, M,levels = [f0],
                 colors=('c',),linestyles=('-',),linewidths=(2,))
plt.clabel(CS, fontsize=9, inline=1)
M_legend0 = mlines.Line2D([], [], color='c', label='$m_k(x_0,y_0)$')
# Solving subproblem
pN = -np.linalg.inv(B0)*g0  #Newton Step
XN = x0 + pN.item(0)
YN = y0 + pN.item(1)
MN = f0 + g0.item(0)*(XN-x0)+g0.item(1)*(YN-y0) + 0.5*B0.item((0,0))*np.multiply((XN-x0),(XN-x0)) + B0.item((0,1))*np.multiply((XN-x0),(YN-y0)) + 0.5*B0.item((1,1))*np.multiply((YN-y0),(YN-y0))
dD = Dmax/nD
xtc = np.zeros(nD)
ytc = np.zeros(nD)
th = np.linspace(0,2*np.pi,nth)
ct = np.cos(th)
st = np.sin(th)
for k in range(1, nD+1):
    Delta = k*dD;
    X = x0 + Delta*ct
    Y = y0 + Delta*st
    M = f0 + g0.item(0)*(X-x0)+g0.item(1)*(Y-y0) + 0.5*B0.item((0,0))*np.multiply((X-x0),(X-x0)) + B0.item((0,1))*np.multiply((X-x0),(Y-y0)) + 0.5*B0.item((1,1))*np.multiply((Y-y0),(Y-y0))
    lowpt = [i for i, v in enumerate(M) if v == min(M)] # minpos = [i for i, v in enumerate(array name here) if v == min(array name here)]
    lowpt = lowpt[0]
    nP = np.linalg.norm(np.subtract([XN, YN],[x0,y0]))
    if (nP>Delta) or (MN > M[lowpt]):
        xtc[k-1] = X[lowpt];
        ytc[k-1] = Y[lowpt];
    else :
        xtc[k-1] = XN;
        ytc[k-1] = YN;
Ax = np.insert(xtc,0,x0)
Ay = np.insert(ytc,0,y0)
plt.plot(Ax, Ay, 'g')
sol_legend = mlines.Line2D([], [], color='g', label='Solution')
plt.legend(handles=[x0y0_legend, f_legend, f_legend0, M_legend, M_legend0, sol_legend],loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
plt.show()
