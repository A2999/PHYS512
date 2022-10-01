import numpy as np
import matplotlib.pyplot as plt

def func(x,y):
    dydx = y/(1+x**2)
    return dydx

def rk4_step(fun,x,y,h):
    k1 = fun(x,y)*h
    k2 = h*fun(x+h/2,y+k1/2)
    k3 = h*fun(x+h/2,y+k2/2)
    k4 = h*fun(x+h,y+k3)
    dy = (k1+2*k2+2*k3+k4)/6
    return y+dy


def rk4_stepd(fun,x,y,h):
    #single step
    k1 = fun(x,y)*h
    k2 = h*fun(x+h/2,y+k1/2)
    k3 = h*fun(x+h/2,y+k2/2)
    k4 = h*fun(x+h,y+k3)
    dy = (k1+2*k2+2*k3+k4)/6
    single = y+dy
    #double
    h2=h/2
    k12 = k1/2 #because of half step
    k22 = h2*fun(x+h2/2,y+k12/2)
    k32 = h2*fun(x+h2/2,y+k22/2)
    k42 = h2*fun(x+h2,y+k32)
    dy2 = (k12+2*k22+2*k32+k42)/6
    y2 = y+dy2
    #calculating the actual double step
    x2=x+h2
    k13 = fun(x2,y2)*h2
    k23 = h2*fun(x2+h2/2,y2+k13/2)
    k33 = h2*fun(x2+h2/2,y2+k23/2)
    k43 = h2*fun(x2+h2,y2+k33)
    dy3 = (k13+2*k23+2*k33+k43)/6
    double = y2+dy3
    
    return double + (double-single)/15

#Comparing the 2 functions

y0=1
steps=200
dsteps=int((4/11)*steps)
xs=np.linspace(-20,20,steps+1)
xd=np.linspace(-20,20,dsteps+1)
hs=(max(xs)-min(xs))/steps
hd=(max(xd)-min(xd))/dsteps
ys=np.zeros(len(xs))
yd=np.zeros(len(xd))
ys[0]=y0
yd[0]=y0

for i in range(len(xs)-1):
    ys[i+1]=rk4_step(func,xs[i],ys[i],hs)
    
for i in range(len(xd)-1):
    yd[i+1]=rk4_stepd(func,xd[i],yd[i],hd)

#Analytical solution
c0 = 1/(np.exp(np.arctan(-20)))
reals = c0*np.exp(np.arctan(xs))
reald = c0*np.exp(np.arctan(xd))

plt.figure(figsize=(8,6))
plt.title('Comparing the accuracy of single step vs double step RK4', fontsize=14)
plt.xlabel('x', fontsize=14)
plt.ylabel('y', fontsize=14)
plt.plot(xd,yd, color='orange', label='RK4 Double Step estimate', linewidth=4)
plt.plot(xs,ys, color='firebrick', label='RK4 estimate')
plt.plot(xs,reals, color='deepskyblue', linestyle='--', label='Real')
plt.legend(fontsize=14, loc='lower right')
plt.show()

errs = np.mean(np.abs(reals-ys))
errd = np.mean(np.abs(reald-yd))
ratio = round(errs/errd,2)

print('Using the double step-size method, we can achieve and accuracy which is {} times better than the single step method.'.format(ratio))

    
    
    