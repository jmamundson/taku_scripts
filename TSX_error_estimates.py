#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 11:49:23 2025

@author: jason
"""

# need to fix so that it pulls the correct resolution from the original images!!!

import rasterio
from osgeo import gdal
import numpy as np
from matplotlib import pyplot as plt
import subprocess
from osgeo import osr

from pyproj import Proj, transform
from pyproj import Transformer

dst_crs = 'EPSG:32608'

import os
import glob
from os.path import basename, isfile


#%%

def importTerraSARX(file):
    '''
    Loads a TerraSAR-X velocity field (geotiff) with rasterio, replaces the minimum data 
    value with NaN, and exports arrays.

    Parameters
    ----------
    file : geotiff filename

    Returns
    -------
    X, Y : coordinates of pixel centers
    data : array of the velocity component
    src_crs : CRS of the data source, used for reprojecting
    '''
    
    with rasterio.open(file) as src:
         data = src.read(1)
         src_crs = src.crs
         
         data[data==-2e9] = np.nan
         height = data.shape[0]
         width = data.shape[1]
         cols, rows = np.meshgrid(np.arange(width), np.arange(height))
         xs, ys = rasterio.transform.xy(src.transform, rows, cols)
         
         # coordinates of pixel centers
         X = np.array(xs)
         Y = np.array(ys)      
     
    return(X, Y, data, src_crs)

#%%

base = '/hdd2/taku/TerraSAR-X/velocities/Release5/Alaska-Taku/'
base = '/hdd2/taku/Sentinel1/Release1-12day/'

dirs = []

for entry in os.scandir(base):
    if entry.is_dir():
        dirs = dirs + [entry.path]
        continue

j = -5

exFile = glob.glob(dirs[j]+'/*_ex_*.tif')[0]
eyFile = glob.glob(dirs[j]+'/*_ey_*.0.tif')[0]
vxFile = glob.glob(dirs[j]+'/*_vx_*.0.tif')[0]
vyFile = glob.glob(dirs[j]+'/*_vy_*.0.tif')[0]


X, Y, ex, src_crs = importTerraSARX(exFile)
_, _, ey, _ = importTerraSARX(eyFile)
_, _, vx, _ = importTerraSARX(vxFile)
_, _, vy, _ = importTerraSARX(vyFile)

vv = np.sqrt(vx**2 + vy**2)

evv = np.sqrt((vx*ex)**2 + (vx*ey)**2)/vv

plt.imshow(evv)
plt.colorbar(label='Error [m a$^{-1}$]')
