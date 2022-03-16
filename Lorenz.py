from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


t_m=100
timepts=10000
x=np.zeros(t_m)
y=np.zeros(t_m)
z=np.zeros(t_m)

sigma, beta, rho = 10, 2.667, 28

x[0], y[0], z[0] = 1.0, 1.0, 1.0
def Lorenz_cont(t,Z,sigma,rho,beta):
    x,y,z=Z
    dx=sigma*(y-x)
    dy=x*(rho-z)-y
    dz=x*y-beta*z
    return dx,dy,dz


L_sol=solve_ivp(Lorenz_cont,(0,t_m),(x[0], y[0], z[0]),args=(sigma,rho,beta),dense_output=True, method='DOP853')

t = np.linspace(0, t_m, timepts)
x,y,z = L_sol.sol(t)




fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_facecolor('k')

s=15
cmap = plt.cm.cool

for i in range(0,timepts-s,s):
    ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=cmap(i/timepts), alpha=0.4)

ax.set_axis_off()
plt.show()
