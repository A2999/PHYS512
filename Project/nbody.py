import numpy as np
import matplotlib.pyplot as plt

def point_charge_V(g, soft_r, b=5):
    x, y = np.meshgrid(np.linspace(-b,b, g), np.linspace(-b,b, g))
    r = np.sqrt(x**2+y**2)
    r[r<soft_r]=soft_r
    V = 1/r
    return V
    
class nbody:
    def __init__(self, n, gridcells, soft_r, x=False, m=False, v=False, plotting=True):
        self.n = n
        if x is False:
            self.x = np.random.randn(n,2)
        else:
            if x.shape[1]!=2:
                print('Position must be a Nx2 array')
                assert(1==0)
            elif int(n)!=len(x):
                print("Number of particles doesn't match number of initial positions.")
                assert(1==0)
            self.x = x
        if m is False:
            self.m = np.ones(n)
        else:
            if int(n)!=len(m):
                print("Number of particles doesn't match number of initial masses.")
                assert(1==0)
            self.m = m
        if v is False:
            self.v = np.zeros([n,2])
        else:
            if v.shape[1]!=2:
                print('Position must be a Nx2 array')
                assert(1==0)
            elif len(x)!=len(v):
                print("Number of initial positions doesn't match number of initial velocities.")
                assert(1==0)
            self.v = v
        self.f = np.empty([self.n,2])
        self.g = gridcells
        self.soft = soft_r
        self.p = np.zeros([self.g,self.g])
            
    def density(self):
        density = np.histogram2d(self.x[:,0],self.x[:,1], bins=self.g, weights=self.m)
        self.p = density[0]
        return self.p
    
    def potential(self):
        density = self.density()
        temp = density.copy()
        temp = np.pad(temp,(0,temp.shape[0]))
        tempft = np.fft.rfft2(temp)
        
        point = point_charge_V(2*self.g, self.soft)
        pointft = np.fft.rfft2(point)
        
        temp = np.fft.irfft2(tempft*pointft)
        
        self.pot = temp[:density.shape[0],:density.shape[1]]
        return self.pot
    
    def force(self):
        V = self.potential()
        ax, ay = np.gradient(V)
        dx = (max(self.x[:,0])-min(self.x[:,0]))/self.g
        dy = (max(self.x[:,1])-min(self.x[:,1]))/self.g
        if dx==0:
            dx=1e-9
        if dy==0:
            dy=1e-9
        for i in range(len(self.x)):
            indx, indy = int((self.x[i,0]-min(self.x[:,0]))//dx), int((self.x[i,1]-min(self.x[:,1]))//dy)
            if indx==self.g: #This has to be done so that the points laying
                indx=indx-1  #on the far edge are included in the previous
            if indy==self.g: #gridcell. Ex: g=5, then c0=[0,1),c1=[1,2),
                indy=indy-1  #c2=[2,3),c3=[3,4),c4=[4,5]
            self.f[i] = -np.array([ax[indx,indy],ay[indx,indy]])*self.m[i]
        return self.f
    
    def leapfrog(self, dt):
        self.x = self.x + self.v*dt
        self.force()
        self.v = self.v + self.f*dt


#parts=nbody(2, 100, 0.001, x=np.array([[0.4,0.5],[0.6,0.5]]))
parts=nbody(100,100,0.001)

plt.ion()
nstep=300
#kk=np.zeros(nstep)
#pp=np.zeros(nstep)
for i in range(nstep):
    plt.clf()
    plt.plot(parts.x[:,0],parts.x[:,1],'.')
    #plt.imshow(parts.p)
    plt.xlim(-3,3)
    plt.ylim(-3,3)
    plt.pause(0.001)
    parts.leapfrog(0.01)
    #kk[i],pp[i]=parts.leapfrog(0.01)