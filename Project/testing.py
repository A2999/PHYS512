# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:57:32 2022

@author: adesr
"""

import numpy as np
import matplotlib.pyplot as plt




print(os.getcwd())




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