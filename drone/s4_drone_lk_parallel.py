#!/usr/bin/env python

import os
import os.path as osp
import subprocess
import glob
import platform

# import pandas as pd
import datetime as dt
import numpy as np

import cv2
from PIL import Image

import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt 
from matplotlib import collections as mc

import multiprocessing as mp

from imports.stop_watch import stop_watch 
import imports.utilities as util


from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import griddata

from osgeo import gdal

from mpl_toolkits.axes_grid1 import make_axes_locatable
#%%

''' 
script to track features using Lucas-Kanade

VARIABLES:
    
source_workspace:
    source workspace, e.g., '/hdd/glaciome/lab_experiments/11_29 run/'
    location of images that will be processed
    
target_workspace:
    workspace for results, e.g., '/hdd/glaciome/lab_experiments/11_29 run/results/'
    location for saving npz files and figures
      


track_len: 
    track length in number of images, e.g. 2. '2' means corners are 
    detected on first image of an image triplet and tracked over two subsequent images: 
    [0,1,2][2,3,4][4,5,6]...
    track_len should be at least 2
    greater track_lens yield increasingly reliable yet thin fields
    greater track_lens work best for photos taken at high rates

startlist:
    list to control at what rate new corners are detected (e.g., every image, every other image)
    in case startlist = [0] and track_len = 2, processing is as follows: [0,1,2][2,3,4]... end => no overlap between velocity fields
    in case startlist = [0, 1] and track_len = 2, processing is a follows: [0,1,2][2,3,4]... end
    ... start over... [1,2,3][3,4,5] ... done
    another example: in case startlist = [0, 2] and track_len = 4, processing is a follows:
    [0,1,2,3,4][4,5,6,7,8]... end ... start over... [2,3,4,5,6][6,7,8,9,10] ... done
    note: numbers in startlist must be integers that are smaller than track_len
    a startlist value of 0 means to start at the first image and loop through all images with no overlap.
    a startlist value of 1 means to start at the second image and loop through all images with no overlap.
    longer startlists provide additional (potentially new) vectors, but also reflect velocities 
    multiple times; also processing will take longer
                        
mask_switch
    switch [1 or 0] to control whether the image is masked
    *mask not yet created, so set this to 0 for the time being
    
plot_switch
    switch [1 or 0] to control whether plots with tracking results are created in the
    subfolders target_workspace/camname/plots
    
movie_switch
    switch [1 or 0] to control whether above plots are animated, 
    one movie per day under target_workspace/camname/plots 
    uses mencoder and works on linux only
    
NOTE: 
    parameters of the Shi-Tomasi corner detector and the Lucas-Kanade tracker 
    are defined in the lucaskanade_tracking function      
'''   
 

def main():
        
    source_workspace = r'/hdd/taku/uav/surveys/03_hillshades/'
    target = r'/hdd/taku/uav/surveys/04_lk/'
            
    track_len = 1
    track_len_sec = 1 # set time difference between photos (this script figures out the time within lucaskanade_tracking)
    
    startlist = [0]
    mask_switch = 0
    
    plot_switch = 1
    movie_switch = 0
    
    #--------------------------------------------------------------------------
    
    
    # start timer
    sw = stop_watch()
    sw.start()
                                         
    lucaskanade_tracking(source_workspace, target, track_len, track_len_sec,
                      startlist, mask_switch, plot_switch, movie_switch)
                            
    # stop timer        
    sw.print_elapsed_time()
    sw.stop() 
    
    
#%%
def lucaskanade_tracking(ws_source, ws_target, track_len, track_len_sec, 
                         startlist, mask_switch, plot_switch, movie_switch):
    
    #-----------------------------------------------
    # parameters for Shi-Tomasi corner detector
    feature_params = dict(maxCorners = 50000000, 
                          qualityLevel = 0.005, 
                          minDistance = 10, 
                          blockSize = 5)    
    
    # parameters for the Lucas-Kanade tracker                       
    lk_params = dict(winSize = (100,100), 
            maxLevel = 2,
            criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 100, 0.03))
   

    # list the imagery in the original folder  
    imagelist = sorted(glob.glob(ws_source + '/*.tif'))
    
    # use following two lines to select specific files from imagelist
    index = [0, 1, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15 ]
    imagelist = [imagelist[j] for j in index]
    
    # avoid folders with image_numbers smaller than track length (would throw errors)
    if len(imagelist) > track_len:
        
        # use gdal to open geotiff files
        ds = gdal.Open(imagelist[0])        
                
        filename = os.path.basename(imagelist[0])
        
        # create time variable for first image
        Y = int(filename[9:13])
        m = int(filename[13:15])
        d = int(filename[15:17])
        H = 00 
        M = 00 
        S = 00
        
        stop_time = dt.datetime(Y,m,d,H,M,S) # this gets renamed after the first time step...
        
                
        mask = np.zeros_like(ds.ReadAsArray())
        
        # here is where you would insert code to mask out the area outside the tank...
        # if mask_switch == 1: 
        #     
        #     #mask[mask1==1] = 255
        # else: # use no mask, set everything to 255
        #     mask[:] = 255
        
        mask[:] = 255 # set to 255 if no mask 
        # empty lists for tracks and track quality 
        # tracks is a list of lists, with each sublist containing the x and y 
        # coordinates of the vertices
        # trackquality contains the pixel distance between tracked and backtracked
        # corners
        tracks = []
        trackquality = []
        
        for start in startlist:        
        
            # loop over photos
            for counter, image in enumerate(imagelist[start:]):
                
                # use gdal to open geotiff files
                ds = gdal.Open(image)        
                frame_gray = ds.ReadAsArray().astype('uint8')
                
                gt = ds.GetGeoTransform()
                ulx = gt[0]
                dx = gt[1]
                uly = gt[3]
                
                lrx = ulx + frame_gray.shape[1]*dx
                lry = uly - frame_gray.shape[0]*dx 
                                
                
                if len(tracks) > 0:
                    
                    filename = os.path.basename(image)
                    print(filename)
                    # this yields the same time stamp as the first image...
                    # create time variable for second image
                    Y = int(filename[9:13])
                    m = int(filename[13:15])
                    d = int(filename[15:17])
                    H = 00 #int(filename_b[13:15])
                    M = 00 #int(filename_b[15:17])
                    S = 00
                    
                    stop_time = dt.datetime(Y,m,d,H,M,S)
            
                    delta_t = stop_time-start_time
            
                    days = delta_t.days + delta_t.seconds/float(86400)
                    
                    # allocate images
                    img0 = prev_gray
                    img1 = frame_gray

                    # high pass filter the hillshades
                    img0 = img0 - gaussian_filter(img0,5,order=0)
                    img1 = img1 - gaussian_filter(img1,5,order=0)
                    
                    # obtain last point from each track
                    p0 = np.float32([tr[-1] for tr in tracks]).reshape(-1, 1, 2)
                    
                    # run tracker
                    p1, st, err = cv2.calcOpticalFlowPyrLK(img0, img1, p0, None, **lk_params)
                    
                    # backtracking for match verification
                    p0r, st, err = cv2.calcOpticalFlowPyrLK(img1, img0, p1, None, **lk_params)
                    
                    # calculate distance between original and backtracked coordinates of p0
                    diff = abs(p0 - p0r).reshape(-1, 2)
                    dist = np.hypot(diff[:, 0], diff[:, 1])
                     
                    # create boolean array for the points with differences < 1 pixel
                    valid = dist < 1
                    
                    new_tracks = []
                    new_trackquality = []
                    
                    # loop over old tracks (tr) and corresponding new extensions (x ,y)
                    # only tracks that can be extented are retained
                    for tr, (x, y), valid_flag, trq, d in zip(tracks, p1.reshape(-1, 2), valid, trackquality, dist):
                        
                        # only add the track extension if flag says 1 (= valid)
                        if valid_flag == 1:
                        
                            # append x, y coordinates to the tracklist
                            tr.append((x, y))
                            trq.append(d)
                            
                            # delete track elements longer than track_len
                            # -1 to account for the fact that tr corresponds to the number of vertices
                            if (len(tr) - 1) > track_len:
                                del tr[0]
                                
                            new_tracks.append(tr)
                            new_trackquality.append(trq)
                    
                    # update the tracklist with new tracks
                    tracks = new_tracks
                    trackquality = new_trackquality 
        
                # save and plot tracking results, determine new corners
                if counter % track_len == 0: # note: "%" is modulus
                                                   
                    if counter > 0 and len(tracks)>0:
                            
                        # calculate velocities and save tracks                                                   
                        tracks = np.array(tracks)
                        
                        # starting point of tracks
                        xi = tracks[:,0,0]
                        yi = tracks[:,0,1]
                                                
                        # ending point of tracks
                        xf = tracks[:,-1,0]
                        yf = tracks[:,-1,1]
                                                
                        # x- and y-components of velocity (pixels/second)
                        u = (xf-xi)/days
                        v = (yf-yi)/days
                        
                        
                        # save tracks
                        npzname = '{}/taku_{}-{}_tracks.npz'.format(ws_target, start_time.strftime('%Y%m%d'), stop_time.strftime('%Y%m%d'))
                        np.savez(npzname, x=xi+ulx, y=uly-yi, u=u, v=-v)#, trackquality = trackquality)
                        
                 
                        # uncomment the following to grid the data without applying any filters
                        # xg = np.arange(0,frame_gray.shape[1])
                        # yg = np.arange(0,frame_gray.shape[0])
                        # X, Y = np.meshgrid(xg, yg)
                        # ug = griddata((xi, yi), u, (X, Y), method='cubic')
                        # vg = griddata((xi, yi), v, (X, Y), method='cubic')
                        # npzname = '{}/{}_{}frames_gridded.npz'.format(ws_target_npz, osp.splitext(osp.basename(image))[0], track_len * track_len_sec)
                        # np.savez(npzname, X=X, Y=Y, ug=ug, vg=vg)
                                
                
                        if plot_switch == 1:
                            
                            # plot_ws = osp.join(ws_target, 'plots') 
                            
                            # if osp.isdir(plot_ws) == 0:
                            #     os.makedirs(plot_ws)                                
                            
                            # turn interactive plotting off
                            plt.ioff() 
                            
                            figsize = (15.0, 15.0 * frame_gray.shape[0] / frame_gray.shape[1])
                            fig, ax = plt.subplots(1, 1, figsize = figsize, facecolor = 'w')
                            
                            img = ax.imshow(frame_gray, cmap = 'gray')
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes('right', size='2%', pad=0.05)                               
                                    
                            # line collection is ideal to print many lines            
                            #ax.add_collection(mc.LineCollection(tracks, color = 'red', alpha = 0.6, linewidth=3))
                            
                            endpoints = np.float32([tr[-1] for tr in tracks]).reshape(-1, 2)
                            
                            #ax.plot(endpoints[:, 0], endpoints[:, 1], '.', color = 'red', ms = 15, alpha = 0.6)
                            speed = np.sqrt(u**2+v**2)
                            
                            v=-v
                            v[v>0] = np.nan
                            v[speed>0.5] = np.nan
                            
                            #ax.quiver(xi,yi,u/speed,-v/speed,speed, scale=100, width=0.001, )
                            h = ax.quiver(xi,yi,u,v,speed, scale=20, width=0.001, clim=[0,0.5] )
                            
                            
                            ax.set_xlim([300, frame_gray.shape[1]-1100])
                            ax.set_ylim([frame_gray.shape[0]-1000, 400])
                                                       
                            ax.set_xticklabels([])
                            ax.set_yticklabels([])
    
    
                            aspect_ratio = (frame_gray.shape[0]-600)/(frame_gray.shape[1]-900)
                            plt.colorbar(h, cax=cax)
                            plt.ylabel('Velocity [m/d]')
                            fig.tight_layout()
                            
                            util.annotatefun(ax, [start_time.strftime('%Y/%m/%d') + '-' + stop_time.strftime('%Y/%m/%d')], 0.03, 0.93, fonts = 22, col = 'white')
                            
                            plotname = '{}/taku_{}-{}_tracks.png'.format(ws_target, start_time.strftime('%Y%m%d'), stop_time.strftime('%Y%m%d'))
                            plt.savefig(plotname, format = 'png', dpi = 80)
                            
                            plt.close(fig)
                        
                            
                            

                    # call Shi-Tomasi corner detector
                    p = cv2.goodFeaturesToTrack(frame_gray, mask = mask, **feature_params)
                                    
                    # reset tracks
                    tracks = []
                    trackquality = []
                    # path string required to save the .npz file with the correct name
                    path = image
                    
                    if p is not None:
                        for x, y in np.float32(p).reshape(-1, 2):
                            tracks.append([(x, y)])
                            trackquality.append([])
                                    
                prev_gray = frame_gray
                start_time = stop_time  
                
        
    
#%%
      
if __name__ == '__main__':
    
    main()


