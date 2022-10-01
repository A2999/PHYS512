import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

#list of all half-lives in seconds
y=60*60*24*365
d=60*60*24
h=60*60
m=60


half_lives = np.array([4.468e9*y, 24.1*d, 6.7*h, 245500*y, 75380*y, 1600*y, 3.8235*d, 3.1*m, 26.8*m, 19.9*m, 164.3e-6, 22.3*y, 5.015*y, 138.376*d, np.inf])

def HL_ODE(x,y,hl=half_lives):
    dydx = np.empty(len(hl))
    #initializing the edge terms
    dydx[0] = -y[0]/hl[0]
    dydx[-1] = y[-2]/hl[-2]
    
    for i in range(len(dydx)-1):
        dydx[i+1] = -y[i+1]/hl[i+1] + y[i]/hl[i]
    return dydx

#initialize the amounts so we start with a pure sample of uranium
y = np.zeros(len(half_lives))
y[0] = 1 

#use Radau because we have a stiff system of equations
ans = integrate.solve_ivp(HL_ODE, [0,half_lives[0]*7], y, t_eval=np.linspace(0, half_lives[0]*7,1000001), method='Radau')
#you can change the t_eval if you plan to run the code,  just have it high to
#provide better resolution for thorium230 vsU234

#Pb206 vs U238
'''
plt.figure(figsize=(8,8))
plt.subplot(2,1,1)
plt.title('Comparing the Proportion of Pb206 and U238 as a function of time', fontsize=14)
plt.plot(ans.t, ans.y[0], label='U238')
plt.plot(ans.t, ans.y[-1], label='Pb206')
plt.ylabel('Proportion of Element', fontsize=14)
plt.legend(fontsize=14)
plt.subplot(2,1,2)
plt.plot(ans.t, ans.y[-1]/ans.y[0])
plt.xlabel('Time (s)', fontsize=14)
plt.ylabel('Ratio of Pb206/U238', fontsize=14)
plt.show()
'''

#Thorium 230 vs U234
plt.figure(figsize=(8,8))
plt.subplot(2,1,1)
plt.title('Comparing the Proportion of Thorium 230 and U234 as a function of time', fontsize=14)
plt.plot(ans.t, ans.y[3], label='U234')
plt.plot(ans.t, ans.y[4], label='Thorium 230')
plt.xlim(-1e13, 0.2e15)
plt.ylabel('Proportion of Element', fontsize=14)
plt.legend(fontsize=14)
plt.subplot(2,1,2)
plt.plot(ans.t, ans.y[4]/ans.y[3])
plt.xlabel('Time (s)', fontsize=14)
plt.ylabel('Ratio of Thorium230/U234', fontsize=14)
plt.xlim(-1e13, 0.2e15)
plt.show()
#'''