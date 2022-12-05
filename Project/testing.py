# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:57:32 2022

@author: adesr
"""

import numpy as np
import matplotlib.pyplot as plt

def point_charge_V(g, soft_r, b=5):
    x, y = np.meshgrid(np.linspace(-b,b, g+1), np.linspace(-b,b, g+1))
    r = np.sqrt(x**2+y**2)
    r[r<soft_r]=soft_r
    V = 1/r
    V=V[:-1,:-1]
    return V
    
class nbody:
    def __init__(self, n, gridcells, soft_r, x=False, m=False, v=False):
        self.n = n
        self.x = np.random.rand(n,2)
        self.m = np.ones(n)
        self.v = np.zeros([n,2])
        self.f = np.zeros([self.n,2])
        self.g = gridcells
        self.soft = soft_r
        self.point = point_charge_V(self.g, self.soft)
        self.pointft = np.fft.rfft2(self.point)           
            
    def density(self):
        density = np.histogram2d(self.x[:,0],self.x[:,1], bins=self.g, range=[[0,1],[0,1]])[0]
        self.p = density
        return self.p
    
    def potential(self):
        density = self.density()
        densityft = np.fft.rfft2(density)
        self.pot = np.fft.fftshift(np.fft.irfft2(densityft*self.pointft))
        
        return self.pot
    
    def force(self):
        V = self.potential()
        fx, fy = np.gradient(V)
        
        dx = 1/self.g
        dy = 1/self.g
        for i in range(len(self.x)):
        
            indx, indy = int(self.x[i,0]/dx), int(self.x[i,1]/dy)
            
            self.f[i] = np.array([fx[indx,indy],fy[indx,indy]])
        return self.f
    
    def leapfrog(self, dt):
        self.x=self.x+dt*self.v
        for i in self.x:
            if i[0]>1:
                i[0]=i[0]-1
            if i[0]<0:
                i[0]=i[0]+1
            if i[1]>1:
                i[1]=i[1]-1
            if i[1]<0:
                i[1]=i[1]+1
        self.v=self.v+self.force()*dt
         

parts=nbody(10000,200,0.05)

"""
rho = parts.density()
pot = parts.potential()
fx, fy = np.gradient(pot)
print(parts.x/(1/parts.g))
plt.imshow(pot)
#"""

"""
plt.subplot(1,2,1)
plt.plot(parts.x[:,0],parts.x[:,1],'.')
plt.subplot(1,2,2)
plt.imshow(parts.point)
plt.colorbar()#"""

#"""
n=100
for i in range(n):
    plt.clf()
    #plt.plot(parts.x[:,0],parts.x[:,1],'.')
    plt.imshow(parts.potential())
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.pause(0.001)
    parts.leapfrog(0.005)
    #"""








































"""
import cv2
import os
import imageio

def plots(particles, nsteps, step='leapfrog', savedir=False):
    if step != 'leapfrog' or step !='rk4':
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
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.leapfrog(0.01)
                plt.savefig('{}lf_{}'.format(savedir,i))
            images = []
            for file in os.listdir(savedir):
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}leapfrog.gif'.format(savedir), images)
    
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
            for i in range(nsteps):
                plt.clf()
                plt.plot(particles.x[:,0],particles.x[:,1],'.')
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.pause(0.001)
                particles.rk4(0.01)
                plt.savefig('{}rk4_{}'.format(savedir,i))
            images = []
            for file in os.listdir(savedir):
                img = cv2.imread(os.path.join(savedir,file))
                images.append(img)
            imageio.mimsave('{}rk4.gif'.format(savedir), images)
"""








"""
#data = np.random.rand(50,2)

def point_charge_V(g, soft_r, b=5):
    x, y = np.meshgrid(np.linspace(-b,b, g), np.linspace(-b,b, g))
    r = np.sqrt(x**2+y**2)
    r[r<soft_r]=soft_r
    V = 1/r
    return V

bins=50

pV = point_charge_V(bins,0.01)
pVft = np.fft.rfft2(pV)
pV2 = point_charge_V(2*bins,0.01)
pV2ft = np.fft.rfft2(pV2)


hist = np.histogram2d(data[:,0],data[:,1], bins=bins)[0]
hist2 = np.zeros([bins*2,bins*2])
hist2[bins//2:3*bins//2, bins//2:3*bins//2] = hist

histft = np.fft.rfft2(hist)
hist2ft = np.fft.rfft2(hist2)

im1 = np.fft.fftshift(np.fft.irfft2(pVft*histft))
im2 = np.fft.fftshift(np.fft.irfft2(pV2ft*hist2ft))[bins//2:3*bins//2, bins//2:3*bins//2]


plt.subplot(1,2,1)
plt.imshow(np.gradient(im1)[0])
plt.subplot(1,2,2)
plt.imshow(np.gradient(im2)[0])"""