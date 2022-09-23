import numpy as np


def integrate_adaptive(fun, a, b, tol, extra=None):
    global adapt_count
    x=np.linspace(a,b,5)
    dx=x[1]-x[0]
    if type(extra)==type(None):
        y=fun(x)
        fx=np.array([[x[0],y[0]],[x[1],y[1]],[x[2],y[2]],[x[3],y[3]],[x[4],y[4]]])
        adapt_count+=5
    else:
        y=np.zeros(5)
        for i in range(len(x)):
            for j in range(len(extra[:,0])):
                if x[i]==extra[j,0]:
                    y[i]=extra[j,1]
                    
        for i in range(len(y)):
            if y[i]==0:
                y[i]=fun(x[i])
                fx=np.vstack((extra,[x[i],y[i]]))
                adapt_count+=1
                
    i3=(y[0]+4*y[2]+y[4])/3*(2*dx)
    i5=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    err=np.abs(i5-i3)
    if err<tol:
        return i5
    else:
        mid=(a+b)/2
        int1=integrate_adaptive(fun,a,mid,tol/2,fx)
        int2=integrate_adaptive(fun,mid,b,tol/2,fx)
        return int1+int2


#Comparing to our inital function
def integrate(fun,a,b,tol):
    global norm_count
    x=np.linspace(a,b,5)
    dx=x[1]-x[0]
    y=fun(x)
    norm_count+=5
    i1=(y[0]+4*y[2]+y[4])/3*(2*dx)
    i2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    myerr=np.abs(i1-i2)
    if myerr<tol:
        return i2
    else:
        mid=(a+b)/2
        int1=integrate(fun,a,mid,tol/2)
        int2=integrate(fun,mid,b,tol/2)
        return int1+int2
    

#Trying with a couple examples
norm_count=0
adapt_count=0    
integrate(np.cos, -5, 5, 1e-6)
integrate_adaptive(np.cos, -5, 5, 1e-6)
cos_diff=norm_count-adapt_count

norm_count=0
adapt_count=0
integrate(np.exp, 0, 10, 1e-6)
integrate_adaptive(np.exp, 0, 10, 1e-6)
exp_diff=norm_count-adapt_count

norm_count=0
adapt_count=0
integrate(np.sin, 0, 10, 1e-6)
integrate_adaptive(np.sin, 0, 10, 1e-6)
sin_diff=norm_count-adapt_count

print('Difference in number of function calls for some simple functions:')
print('cos(x):{}'.format(cos_diff))
print('sin(x):{}'.format(sin_diff))
print('exp(x):{}'.format(exp_diff))