#!/usr/bin/env python

import os
import os.path as osp
import multiprocessing
import subprocess
import glob
import warnings

import pandas as pd
import datetime as dt
import numpy as np


from PIL import Image
from PIL import ImageFile

import matplotlib
#matplotlib.use('Agg')

import matplotlib as mpl
import matplotlib.pyplot as plt 
import matplotlib.path as mplPath 
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import cmocean

from imports.stop_watch import stop_watch 
import imports.utilities as util
import imports.tracking_misc as trm

try:
    warnings.filterwarnings('ignore')
    mpl.rc("font", **{"sans-serif": ["Arial"]})
except:
    pass

from osgeo import gdal, osr

daysYear = 365.25

''' 
script to combine results from multiple cameras, to average velocity fields 
spatially and temporally, and to plot the results

VARIABLES:
   
source_path_head:
    path head portion to the directory with the .npz input files (projected velocity fields)
    files will be searched under source_path_head/camname/source_path_tail
    
camnames: 
    camera(s) to be processed
    single name or list of names, e.g. 'UAS7' or ['UAS7', 'UAS8'] 
    camera name must match a folder name in the 'source_path_head' directory
    it also has to match the camera name in the parameter excel file 'paramfile'
    files will be searched under source_path_head/camname/source_path_tail
    
source_path_tail:
    path tail portion to the directory with the .npz input files (projected velocity fields)
    files will be searched under source_path_head/camname/source_path_tail
    
source_path_photos:
    path to directory with original .jpg photos (used for plotting)
    e.g. '/hdd2/leconte_2018/cameras'
    photos will be searched in daily folders
    source_path_photos/camname/%Y%m%d
    
target_path:
    path to directory that will contain the integrated velocity fields (.npz files)
    and plots
      
paramfile:
    excel file featuring processing parameters, e.g., from camera calibration
    must be placed under ./ (in the same directory as this script)
    e.g., 'parameter_file_2018.xlsx'
    
clockdriftfile:
    excel file featuring clockdrift related info for each camera
    must be placed under ./data/

fjord_outline:
    .npz file with fjord outline coordinates 
    required to create grid across fjord and for plotting
    must be placed under ./data/
    must contain the three arrays 'x, 'y', and 'id', 'id' must be set to 1
    
min_date: 
    start date of period of interest, e.g. 20180904 
max_date: 
    end date of period of interest, e.g. 20180905
    
time_window:
    time window (in hours) for temporal aggregation of trajectories, e.g. 1.0
    24 aggregates over full day
grid_size:
    grid size (in m) used for spatial aggregation of trajectories, e.g. 50
observation_threshold:
    minimum number of trajectories per grid cell, e.g. 10
    
plot_switch:
    switch [0, 1, 2] to control plotting
    0 = create no plots (i.e., create .npz files only)
    1 = create plot with one map showing aggregated velocities 
    2 = create plot with two maps, showing original and aggregated velocities        
        (avoid '2' for time windows > 30 minutes as plotting of original velocities 
        may use too much memory) 
    3 = same as 2, but map entire fjord rather than terminus area only 

speedthreshold_cbar:
    threshold (m/s) limiting upper speed limit of colorbar
    
movie_switch:
    switch [1 or 0] to control whether plots are animated into a movie

n_proc: 
    number of processors used during processing, e.g. 8 

NOTE: 
    adapt plotting variables in the two plotting functions directly    
    
''' 

# clean up script so that there you only have to specify the base directory in one location!
  
def main():
    
    # parameters
    source_path_head = '/hdd/taku/uav/surveys/04_lk/'
    target_path = '/hdd/taku/uav/surveys/05_lk_grids/' 
      
    grid_size = 50
    observation_threshold = 3
    speedthreshold_cbar = 0.5
        
    plot_switch = 1
    
    #--------------------------------------------------------------------------
    if osp.isdir(target_path) == 0:
        os.makedirs(target_path)
    
        
    arguments = (source_path_head, target_path, grid_size, speedthreshold_cbar, observation_threshold, 
                  plot_switch)
    
    utm_to_gridded_utm(arguments)
            
    
    
#%%
        
def savegeotiff(x, y, data, rows, cols, grid_size, geotiff_name):
    
    x = x.reshape((cols,rows)).T
    y = y.reshape((cols,rows)).T
    data = data.reshape((cols,rows)).T
                
    geotransform = (x[0,0], grid_size, 0, y[0,0], 0, -grid_size) 
              
    
    dst_ds = gdal.GetDriverByName('GTiff').Create(geotiff_name, cols, rows, 1, gdal.GDT_Float32)

    dst_ds.SetGeoTransform(geotransform)    # specify coords
    srs = osr.SpatialReference()            # establish encoding
    srs.ImportFromEPSG(32608)                # WGS84 UTM Zone 8N
    dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
    dst_ds.GetRasterBand(1).WriteArray(data*daysYear)   # write r-band to the raster
    dst_ds.FlushCache()                     # write to disk
    dst_ds = None        


def utm_to_gridded_utm(arguments):
    
    (source_path_head, target_path, grid_size, speedthreshold_cbar, observation_threshold, plot_switch) = arguments
    
    # determine directory of current script (does not work in interactive mode)
    # file_path = osp.dirname(osp.realpath(__file__))
    
    #dem_file = 'taku_DEM_20211017_hillshade.tif'

    # create coarse grid across the fjord    
    

    tracks_files = sorted(glob.glob(source_path_head + '*tracks.npz'))
    dem_files = sorted(glob.glob('/hdd/taku/uav/surveys/03_hillshades/*.tif'))
    
    index = [0, 1, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16 ]
    dem_files = [dem_files[j] for j in index]
    
    [polygons_coarse, centerpoints_coarse, indices, topleft, rows, cols] = trm.create_grid(dem_files[0], grid_size)

    
    for j in np.arange(0,len(tracks_files)):    
        tracks = np.load(tracks_files[j])
    
        x = tracks['x']
        y = tracks['y']
        u = tracks['u']
        v = tracks['v']
   
        # adapt velocity shape
        points = np.vstack((x, y)).T
        point_displacements = np.vstack((u, v)).T
    
        # prep empty lists to store results for each grid cell 
        grid_id = []
        coarsegrid_i = []
        coarsegrid_j = []        
        coarsegrid_x = []
        coarsegrid_y = []
        coarsegrid_u = []
        coarsegrid_v = []
        coarsegrid_speed = []
        count = []
        
        polygons_measured = []
        polygons_not_measured = []
        
        # loop over grid and detect velocities within each grid cell    
        for counter, (poly, centerpoint, index) in enumerate(zip(polygons_coarse, centerpoints_coarse, indices)): 
            
            # check which velocities are contained in grid cell
            grid = mplPath.Path(poly).contains_points(points)
            points_selected = points[grid==1]
            displacements_selected = point_displacements[grid==1] 
            
            nr_observations = len(points_selected) 
            
            coarsegrid_i.append(index[0])
            coarsegrid_j.append(index[1])
            coarsegrid_x.append(centerpoint[0]) 
            coarsegrid_y.append(centerpoint[1]) 
        
            if nr_observations > observation_threshold:
                                
                # divide by # of observations to obtain average displacements in m/s
                mean_u = np.sum(displacements_selected[:, 0]) / nr_observations
                mean_v = np.sum(displacements_selected[:, 1]) / nr_observations
                
                coarsegrid_u.append(mean_u)
                coarsegrid_v.append(mean_v) 
                coarsegrid_speed.append(np.hypot(mean_u, mean_v))
                
                count.append(nr_observations)
                grid_id.append(counter)
                
                polygons_measured.append(poly)
                
            else:
                polygons_not_measured.append(poly)
                coarsegrid_u.append(np.nan)
                coarsegrid_v.append(np.nan)

                   
        
        ds = gdal.Open(dem_files[j+1])
        image = ds.ReadAsArray().astype('uint8')
        gt = ds.GetGeoTransform()
        ulx = gt[0]
        dx = gt[1]
        uly = gt[3]
        
        lrx = ulx + image.shape[1]*dx
        lry = uly - image.shape[0]*dx
    
        
        x = np.array(coarsegrid_x)
        y = np.array(coarsegrid_y)
        u = np.array(coarsegrid_u)
        v = np.array(coarsegrid_v)
        speed = np.sqrt(u**2+v**2)
        
        
        npzname = '{}/{}_gridded.npz'.format('/hdd/taku/uav/surveys/05_lk_grids', osp.splitext(osp.basename(tracks_files[j]))[0]) 
        np.savez(npzname, x=x, y=y, u=u, v=v)
       
        
        #-------- plot figure ---------
        figsize = (15.0, 15.0 * image.shape[0] / image.shape[1])
        fig, ax = plt.subplots(1, 1, figsize = figsize, facecolor = 'w')
        
        img = ax.imshow(image, cmap = 'gray', extent=(ulx,lrx,lry,uly))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.05)                               
        
        h = ax.quiver(x,y,u,v,speed, scale=10, width=0.001, clim=[0,0.5] )
        
        
        ax.set_xlim([ulx, lrx-1100])
        #ax.set_ylim([lry+400, uly-1000])
        ax.set_ylim([lry+1000, uly-400])
    
                               
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    
    
        aspect_ratio = (image.shape[0]-600)/(image.shape[1]-900)
        plt.colorbar(h, cax=cax)
        plt.ylabel('Velocity [m/d]')
        fig.tight_layout()
            
                  
        plotname = '{}/{}_gridded.png'.format('/hdd/taku/uav/surveys/05_lk_grids', osp.splitext(osp.basename(tracks_files[j]))[0]) 
        plt.savefig(plotname, format = 'png', dpi = 80)
        
        plt.close()                    
    
    

            
        geotiff_name = '{}/{}_{}.tif'.format('/hdd/taku/uav/surveys/05_lk_grids', osp.basename(tracks_files[j]).rstrip('_tracks.npz'), 'vx')
        savegeotiff(x, y, u, rows, cols, grid_size, geotiff_name)        
        
        geotiff_name = '{}/{}_{}.tif'.format('/hdd/taku/uav/surveys/05_lk_grids', osp.basename(tracks_files[j]).rstrip('_tracks.npz'), 'vy')
        savegeotiff(x, y, v, rows, cols, grid_size, geotiff_name)
        
        geotiff_name = '{}/{}_{}.tif'.format('/hdd/taku/uav/surveys/05_lk_grids', osp.basename(tracks_files[j]).rstrip('_tracks.npz'), 'vv')
        savegeotiff(x, y, np.sqrt(u**2+v**2), rows, cols, grid_size, geotiff_name)
        
    


    
#%%
    
if __name__ == '__main__':
    
    main()
