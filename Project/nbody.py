import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import imageio

def point_charge_V(g, soft_r, b=5):
    x, y = np.meshgrid(np.linspace(-b,b, g), np.linspace(-b,b, g))
    r = np.sqrt(x**2+y**2)
    r[r<soft_r]=soft_r
    V = 1/r
    return V

def get_derivs(xv, n, fun):
    x = xv[:n,:]
    v = xv[n:,:]
    f = fun(x)
    return np.vstack([v,f])
    
class nbody:
    def __init__(self, n, gridcells, soft_r, x=False, m=False, v=False, periodic=True):
        self.n = n
        if x is False:
            self.x = np.random.rand(n,2)
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
        self.f = np.zeros([self.n,2])
        self.g = gridcells
        self.soft = soft_r
        self.periodic = periodic
            
    def density(self, rk4=False):
        if rk4 is False:
            density = np.histogram2d(self.x[:,0],self.x[:,1], bins=self.g, weights=self.m)[0]
        else:
            density = np.histogram2d(rk4[:,0],rk4[:,1], bins=self.g, weights=self.m)[0]
        
        if self.periodic is True:
            self.p = density
        else:
            self.p = np.zeros([self.g*2,self.g*2])
            self.p[self.g//2:3*self.g//2, self.g//2:3*self.g//2] = density
        return self.p
    
    def potential(self, rk4=False):
        if rk4 is False:
            density = self.density()
        else:
            density = self.density(rk4)
        
        densityft = np.fft.rfft2(density)
        
        if self.periodic is True:
            point = point_charge_V(self.g, self.soft)            
            pointft = np.fft.rfft2(point)
            self.pot = np.fft.fftshift(np.fft.irfft2(densityft*pointft))
        else:
            point = point_charge_V(2*self.g, self.soft)            
            pointft = np.fft.rfft2(point)
            big = np.fft.fftshift(np.fft.irfft2(densityft*pointft))
            self.pot = big[self.g//2:3*self.g//2, self.g//2:3*self.g//2]
            
        return self.pot
    
    def force(self, rk4=False):
        if rk4 is False:
            V = self.potential()
        else:
            V = self.potential(rk4)
        fx, fy = np.gradient(V)
        dx = (max(self.x[:,0])-min(self.x[:,0]))/self.g
        dy = (max(self.x[:,1])-min(self.x[:,1]))/self.g
        if dx==0:
            dx=1e-9
        if dy==0:
            dy=1e-9
        for i in range(len(self.x)):
            indx, indy = int((self.x[i,0]-min(self.x[:,0]))//dx), int((self.x[i,1]-min(self.x[:,1]))//dy)
            #print(indx, indy)
            if indx==self.g: #This has to be done so that the points laying
                indx=indx-1  #on the far edge are included in the previous
            if indy==self.g: #gridcell. Ex: g=5, then c0=[0,1),c1=[1,2),
                indy=indy-1  #c2=[2,3),c3=[3,4),c4=[4,5]
            self.f[i] = -np.array([fx[indx,indy],fy[indx,indy]])
        return self.f
    
    def leapfrog(self, dt):
        self.x=self.x+dt*self.v
        self.v=self.v+self.force()*dt
        
    def rk4(self, dt):
        xv = np.vstack([self.x,self.v])
        k1 = get_derivs(xv, self.n, self.force)
        k2 = get_derivs(xv+k1*dt/2, self.n, self.force)
        k3 = get_derivs(xv+k2*dt/2, self.n, self.force)
        k4 = get_derivs(xv+k3*dt, self.n, self.force)
        
        tot=(k1+2*k2+2*k3+k4)/6
        
        self.x = self.x + tot[:self.n,:]*dt
        self.v = self.v + tot[self.n:,:]*dt


def plots(particles, nsteps, step='leapfrog', savedir=False):
    if step != 'leapfrog' and step !='rk4':
        print('Only permitted step mehtods are "leapfrog" and "rk4"')
        assert(1==0)
    if step =='leapfrog':
        if savedir is False:
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.leapfrog(0.01)
        else:
            filenames = []
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.leapfrog(0.01)
                plt.savefig('{}/lf_{}.png'.format(savedir,i))
                filenames.append('lf_{}.png'.format(i))
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/leapfrog.gif'.format(savedir), images)
    
    if step =='rk4':
        if savedir is False:
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.rk4(0.01)
        else:
            filenames = []
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.rk4(0.01)
                plt.savefig('{}/rk4_{}.png'.format(savedir,i))
                filenames.append('rk4_{}.png'.format(i))
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/rk4.gif'.format(savedir), images)
         

parts=nbody(100,100,0.001, periodic=True)

plots(parts, 50, step='leapfrog', savedir='images/leapfrog_periodic')