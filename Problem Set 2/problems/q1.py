import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#Making our integration function
def integrate(fun,a,b,tol):
    #making our 5 points and dx
    x=np.linspace(a,b,5)
    dx=x[1]-x[0]
    y=fun(x)
    
    #do the 3-point integral
    i1=(y[0]+4*y[2]+y[4])/3*(2*dx)
    #do the 5-point integral
    i2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    #Finding the error between the 3 & 5 point integral
    myerr=np.abs(i1-i2)
    #if the error is smaller than our tolerance, then we're set
    if myerr<tol:
        return i2
    #Otherwise we need to split the section and try again
    else:
        mid=(a+b)/2
        int1=integrate(fun,a,mid,tol/2)
        int2=integrate(fun,mid,b,tol/2)
        return int1+int2
    
#Setting up the integral for the electric field of the spherical shell
def func(theta):
    return (sigma*R**2/(2*eps_0))*(z-R*np.cos(theta))*np.sin(theta)/(R**2+z**2-2*R*z*np.cos(theta))**(3/2)

#defining all the variables in our function
R=1
sigma=1
eps_0=1

x=np.linspace(0,2*np.pi,2001)
y=np.zeros(len(x))
y2=np.zeros(len(x))

#Integrating for multiple values of z
for i in range(len(x)):
    z=x[i]
    y[i]=integrate(func,0,np.pi,1e-6)
    y2[i]=quad(func, 0, np.pi)[0]

plt.figure(figsize=(8,6))
plt.title('The Electric Field as a function of distance of a spherical shell of radius {}'.format(R))
plt.plot(x,y, label='integrate function', color='magenta')
plt.plot(x,y2, label='scipy.quad', linestyle='--', color='mediumseagreen')
plt.legend(fontsize=12)
plt.show()

