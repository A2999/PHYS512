import numpy as np


x=np.linspace(-1,1,1001)
#transposing the range x=[0.5,1] to x=[-1,1]
y=np.log2((x+3)/4)

deg=7

coeffs=np.polynomial.chebyshev.chebfit(x,y, deg)

chebs=np.empty([len(x),deg+1])

chebs[:,0]=1
chebs[:,1]=x
for i in range(1,deg):
    chebs[:,i+1]=2*x*chebs[:,i]-chebs[:,i-1]
    
yt = chebs@coeffs

plt.plot(x,y-yt)
print(max(abs(y-yt)))

ln=yt/np.log2(np.e)


def mylog2(x, coeffs, deg):
    num=np.frexp(x)
    a=num[0]*4-3
    
    chebs1=np.empty(deg+1)
    chebs1[0]=1
    chebs1[1]=a
    for i in range(1,deg):
        chebs1[i+1]=2*a*chebs1[i]-chebs1[i-1] 
    y = chebs1@coeffs+num[1]
    
    num2=np.frexp(np.e)
    e=num2[0]*4-3
    chebs2=np.empty(deg+1)
    chebs2[0]=1
    chebs2[1]=e
    for i in range(1,deg):
        chebs2[i+1]=2*e*chebs2[i]-chebs2[i-1] 
    e2 = chebs2@coeffs+num2[1]
    
    return y/e2

print(mylog2(4, coeffs, deg))
print(np.log(4))