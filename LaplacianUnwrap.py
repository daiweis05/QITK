# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 18:01:03 2015

@author: liw5
"""
import numpy as np

def uwnrap_LP(array,voxelsize):
    ArrayDim=array.shape
    FOV=np.array(voxelsize[0:3])*np.array(ArrayDim[0:3])
    FOV=FOV/np.max(FOV)
    
    nx=np.arange(-ArrayDim[1]/2, ArrayDim[1]/2, 1)/FOV[0]
    ny=np.arange(-ArrayDim[0]/2, ArrayDim[0]/2, 1)/FOV[1]
    nz=np.arange(-ArrayDim[2]/2, ArrayDim[2]/2, 1)/FOV[2]
    xx, yy, zz = np.meshgrid(nx, ny,nz)
    k2=xx**2+yy**2+zz**2
    k2= np.fft.fftshift(k2)
    array1=np.cos(array)*np.fft.ifftn(k2*np.fft.fftn(np.sin(array)))
    array1-=np.sin(array)*np.fft.ifftn(k2*np.fft.fftn(np.cos(array)))
    array1=np.fft.fftn(array1)
    k2[0,0,0]=1
    array1=array1/k2
    array1[0,0,0]=1
    array1=np.real(np.fft.ifftn(array1))
    finished=1    
    yield array1
    yield finished

