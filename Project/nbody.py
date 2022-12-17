import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import imageio


def point_charge_V(g, soft_r, b=5):
    x, y = np.meshgrid(np.linspace(-b,b, g+1), np.linspace(-b,b, g+1))
    r = np.sqrt(x**2+y**2)
    r[r<soft_r]=soft_r
    V = 1/r
    V=V[:-1,:-1]
    return V

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
        self.mstack = np.vstack((self.m,self.m)).T
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
        self.g = gridcells
        self.soft = soft_r
        self.periodic = periodic
        self.pot = np.zeros([self.g,self.g])
        if periodic is True:
            self.point = point_charge_V(self.g, self.soft)
            self.pointft = np.fft.rfft2(self.point)
        else:
            self.point = point_charge_V(2*self.g, self.soft)            
            self.pointft = np.fft.rfft2(self.point)
            
            self.lostK = 0
            self.lostU = 0
            
    def density(self, rk4=False):
        if rk4 is False:
            density = np.histogram2d(self.x[:,0],self.x[:,1], bins=self.g, weights=self.m, range=[[0,1],[0,1]])[0]
        else:
            density = np.histogram2d(rk4[:,0],rk4[:,1], bins=self.g, weights=self.m, range=[[0,1],[0,1]])[0]

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
            self.pot = np.fft.fftshift(np.fft.irfft2(densityft*self.pointft))
        else:
            big = np.fft.fftshift(np.fft.irfft2(densityft*self.pointft))
            self.pot = big[self.g//2:3*self.g//2, self.g//2:3*self.g//2]
            
        return self.pot
    
    def force(self, rk4=False):
        if rk4 is False:
            V = self.potential()
        else:
            V = self.potential(rk4)
        fx, fy = np.gradient(V)
        dx = 1/self.g
        dy = 1/self.g
        self.f = np.empty([len(self.x),2])
        for i in range(len(self.x)):
            indx, indy = int(self.x[i,0]/dx), int(self.x[i,1]/dy)
            if indx == self.g:
                indx=indx-1
            if indy == self.g:
                indy=indy-1
            self.f[i] = np.array([fx[indx,indy],fy[indx,indy]])
        return self.f
    
    def leapfrog(self, dt):
        self.x=self.x+dt*self.v
        if self.periodic is True:
            for i in self.x:
                if i[0]>1:
                    i[0]=i[0]-1
                if i[0]<0:
                    i[0]=i[0]+1
                if i[1]>1:
                    i[1]=i[1]-1
                if i[1]<0:
                    i[1]=i[1]+1
        else:
            problems = []
            for i in range(len(self.x)):
                if self.x[i,0]>1 or self.x[i,1]>1 or self.x[i,0]<0 or self.x[i,1]<0:
                    problems.append(i)
            if len(problems)!=0:
                for i in problems:
                    self.lostK = self.lostK + np.sum(0.5*self.mstack[i]*self.v[i]**2)
                    self.lostU = self.lostU + np.sum(self.pot)/self.n
                self.x = np.delete(self.x, problems, 0)
                self.v = np.delete(self.v, problems, 0)
                self.m = np.delete(self.m, problems)
                self.mstack = np.delete(self.mstack, problems, 0)
        self.v=self.v+self.force()*dt
        
    def get_derivs(self,xx):
        nn=xx.shape[0]//2
        x=xx[:nn,:]
        v=xx[nn:,:]
        f=self.force(x)
        return np.vstack([v,f])
    
    def rk4(self, dt):
        print(np.shape(self.x), np.shape(self.v))
        xv = np.vstack([self.x,self.v])
        k1 = self.get_derivs(xv)
        k2 = self.get_derivs(xv+k1*dt/2)
        k3 = self.get_derivs(xv+k2*dt/2)
        k4 = self.get_derivs(xv+k3*dt)
        
        tot=(k1+2*k2+2*k3+k4)/6
        v = tot[:len(self.x),:]
        f=tot[len(self.x):,:]
        
        self.x = self.x + v*dt
        if self.periodic is True:
            for i in self.x:
                if i[0]>1:
                    i[0]=i[0]-1
                if i[0]<0:
                    i[0]=i[0]+1
                if i[1]>1:
                    i[1]=i[1]-1
                if i[1]<0:
                    i[1]=i[1]+1
        else:
            problems = []
            for i in range(len(self.x)):
                if self.x[i,0]>1 or self.x[i,1]>1 or self.x[i,0]<0 or self.x[i,1]<0:
                    problems.append(i)
            if len(problems)!=0:
                for i in problems:
                    self.lostK = self.lostK + np.sum(0.5*self.mstack[i]*self.v[i]**2)
                    self.lostU = self.lostU + np.sum(self.pot)/self.n
                self.x = np.delete(self.x, problems, 0)
                self.v = np.delete(self.v, problems, 0)
                self.m = np.delete(self.m, problems)
                self.mstack = np.delete(self.mstack, problems, 0)
                f = np.delete(f, problems, 0)
        self.v = self.v + f*dt


def plots(particles, nsteps, dt, step='leapfrog', savedir=False, Vplots=False, trackE=False):
    if step != 'leapfrog' and step !='rk4':
        print('Only permitted step mehtods are "leapfrog" and "rk4"')
        assert(1==0)
    plt.figure(figsize=(6,6))
    if savedir is not False:
        filenames = []
    if trackE is not False:
        K = np.zeros(nsteps)
        U = np.zeros(nsteps)
    if step =='leapfrog' and Vplots is not True:
        for i in range(nsteps):
            if trackE is not False:
                K[i] = np.sum(0.5*particles.mstack*particles.v**2)+particles.lostK
                U[i] = np.sum(particles.potential(particles.x))+particles.lostU
            plt.clf()
            plt.plot(particles.x[:,0],particles.x[:,1],'.')
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.pause(0.001)
            particles.leapfrog(dt)
            if savedir is not False:
                plt.savefig('{}/lf_{}.png'.format(savedir,i))
                filenames.append('lf_{}.png'.format(i))
        if savedir is not False:
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/leapfrog.gif'.format(savedir), images)
        
    if step =='leapfrog' and Vplots is True:
        for i in range(nsteps):
            if trackE is not False:
                K[i] = np.sum(0.5*particles.mstack*particles.v**2)+particles.lostK
                U[i] = np.sum(particles.potential(particles.x))+particles.lostU
            plt.clf()
            plt.pcolormesh(particles.pot.T)
            plt.pause(0.001)
            particles.leapfrog(dt)
            if savedir is not False:
                plt.savefig('{}/lf_V_{}.png'.format(savedir,i))
                filenames.append('lf_V_{}.png'.format(i))
        if savedir is not False:
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/leapfrog_V.gif'.format(savedir), images)
    
    if step =='rk4' and Vplots is not True:
        for i in range(nsteps):
            if trackE is not False:
                K[i] = np.sum(0.5*particles.mstack*particles.v**2)+particles.lostK
                U[i] = np.sum(particles.potential(particles.x))+particles.lostU
            plt.clf()
            plt.plot(particles.x[:,0],particles.x[:,1],'.')
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.pause(0.001)
            particles.rk4(dt)
            if savedir is not False:
                plt.savefig('{}/rk4_{}.png'.format(savedir,i))
                filenames.append('rk4_{}.png'.format(i))
        if savedir is not False:
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/rk4.gif'.format(savedir), images)
    
    if step =='rk4' and Vplots is True:
        for i in range(nsteps):
            if trackE is not False:
                K[i] = np.sum(0.5*particles.mstack*particles.v**2)+particles.lostK
                U[i] = np.sum(particles.potential(particles.x))+particles.lostU
            plt.clf()
            plt.pcolormesh(particles.potential(particles.x).T)
            plt.pause(0.001)
            particles.rk4(dt)
            if savedir is not False:
                plt.savefig('{}/rk4_V_{}.png'.format(savedir,i))
                filenames.append('rk4_V_{}.png'.format(i))
        if savedir is not False:
            images = []
            for file in filenames:
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}/rk4_V.gif'.format(savedir), images)
    if trackE is not False:
        return K, U, K+U
         
circlex = np.array([[0.4,0.5],[0.6,0.5]])
circlev = np.array([[0,0.05],[0,-0.05]])
nparts = 2000
xs = np.random.rand(nparts,2)/4+0.375

parts=nbody(nparts,100,0.5, x=xs, periodic=False)
K, U, E = plots(parts, 100, 0.01, step='rk4', Vplots=False, trackE=True)