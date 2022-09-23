import numpy as np


def integrate_adaptive(fun, a, b, tol, extra=None):
    #determining the x coords of our 5 points and dx
    x=np.linspace(a,b,5)
    dx=x[1]-x[0]
    
    #Checking if we have an argument for 'extra' to
    #determine whether this is the first time going
    #through the function
    if type(extra)==type(None):
        #Since this is the first time, we need to
        #calculate all y values
        y=fun(x)
        #save the x & y values in an array
        fx=np.array([[x[0],y[0]],[x[1],y[1]],[x[2],y[2]],[x[3],y[3]],[x[4],y[4]]])
    else:
        #if this isn't the first time we go through
        #the function, we first start by making an
        #'empty' array to which we will save our
        #function values
        y=np.zeros(5)
        #Here we go through the x values we want to
        #find f(x) for and check extra to see if we
        #have already calculated f(x) for said x
        for i in range(len(x)):
            for j in range(len(extra[:,0])):
                #If we have already, then no need to
                #calculate again, simply pull the
                #f(x) value from our array
                if x[i]==extra[j,0]:
                    y[i]=extra[j,1]
        #If we haven't then calculate f(x)
        #and add it to our array
        for i in range(len(y)):
            if y[i]==0:
                y[i]=fun(x[i])
                fx=np.vstack((extra,[x[i],y[i]]))
    
    #Then its like the original integrate function               
    
    #3-point integral
    i3=(y[0]+4*y[2]+y[4])/3*(2*dx)
    #5-point integral
    i5=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/3*dx
    
    err=np.abs(i5-i3)
    if err<tol:
        return i5
    else:
        mid=(a+b)/2
        int1=integrate_adaptive(fun,a,mid,tol/2,fx)
        int2=integrate_adaptive(fun,mid,b,tol/2,fx)
        return int1+int2