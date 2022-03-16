from scipy.integrate import solve_ivp
import numpy as np
t_m=100
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

def Lorenz_discr(x,y,z,n,sigma,rho,beta):
    for i in n:
        x[i+1]=x[i]+schrittweite*(sigma*(y[i]-x[i]))
        y[i+1]=y[i]+schrittweite*(x[i]*(rho-z[i])-y[i])
        z[i+1]=z[i]+schrittweite*(x[i]*y[i]-beta*z[i])
    return x,y,z

print('info', sigma,rho,beta)
L_sol=solve_ivp(Lorenz_cont,(0,t_m),(x[0], y[0], z[0]),args=(sigma,rho,beta),dense_output=True, method='DOP853')
