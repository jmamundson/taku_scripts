#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:57:30 2025

@author: jason
"""

import numpy as np
from matplotlib import pyplot as plt

import rasterio
import glob
import datetime as dt

import os
import subprocess

#%%
year = '2023'

base = '/hdd/taku/uav/surveys/'

#%% load gridded velocities
vxFiles = sorted(glob.glob(base + '05_lk_grids/*vx.tif'))
vyFiles = sorted(glob.glob(base + '05_lk_grids/*vy.tif'))

if year=='2023':
    index = [5, 6]
    vxFiles = [vxFiles[x] for x in index]
    vyFiles = [vyFiles[x] for x in index]

# load DEMs into one stack
vx_dict = {}
vy_dict = {}

for j in np.arange(0,len(vxFiles)):
    with rasterio.open(vxFiles[j]) as src:
        vx_dict[j] = src.read(1)
        transform = src.get_transform()
        bounds = src.bounds    
    with rasterio.open(vyFiles[j]) as src:
        vy_dict[j] = src.read(1)
        transform = src.get_transform()
        bounds = src.bounds    
        
vx_keys = vx_dict.keys()
vx = np.stack([vx_dict[x] for x in vx_keys], axis=0)

vy_keys = vy_dict.keys()
vy = np.stack([vy_dict[x] for x in vy_keys], axis=0)


# use this to resample the DEM to the same resolution
dx = transform[1]
dy = transform[-1]


#%% load DEMS
# need to be careful with sign when taking derivatives!

dems = sorted(glob.glob(base + '02_takuGrid/*DEM*tif'))

if year=='2023':
    index = [6, 10, 11]
    dems = [dems[x] for x in index]

filenames = [os.path.basename(x) for x in dems]

# resample DEMs so they have the same grid as the velocity arrays
for j in np.arange(0,len(dems)):
    subprocess.run(['gdal_translate', '-tr', str(dx), str(dy), '-r', 'cubicspline', dems[j], filenames[j]])


dates = [dt.datetime.strptime(x[9:17],'%Y%m%d') for x in filenames]
deltaT = np.diff(dates)
deltaT = np.array([x.days for x in deltaT]) # time difference in days

# load DEMs into one stack
dems_dict = {}

for j in np.arange(0,len(filenames)):
    with rasterio.open(filenames[j]) as src:
        dems_dict[j] = src.read(1)
        transform = src.get_transform()
        bounds = src.bounds


dems_keys = dems_dict.keys()
dems = np.stack([dems_dict[x] for x in dems_keys], axis=0)
dems[dems<0] = np.nan

dh = np.diff(dems, axis=0)
dhdt = np.array([dh[x,:,:]/deltaT[x] for x in np.arange(0,2)])

dhdx = np.gradient(dems, dx, axis=2)
dhdy = np.gradient(dems, -dy, axis=1)

# compute mean slope during each time interval
dhdx = (dhdx[:-1,:,:]+dhdx[1:,:,:])/2
dhdy = (dhdy[:-1,:,:]+dhdy[1:,:,:])/2

#%% dh/dt + vx*dh/dx+vy*dh/dy

vDownslope = vx*dhdx + vy*dhdy

plt.subplot(1,2,1)
plt.imshow(dhdt[0,:,:]+vDownslope[0,:,:], vmin=-20, vmax=20)

plt.subplot(1,2,2)
plt.imshow(dhdt[1,:,:]+vDownslope[1,:,:], vmin=-20, vmax=20)

        
        
