from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

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

plt.plot(t, z.T)
plt.show()
