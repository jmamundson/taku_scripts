import rasterio
import gdal
import numpy as np
from matplotlib import pyplot as plt
import subprocess
from osgeo import osr

from pyproj import Proj, transform
from pyproj import Transformer

dst_crs = 'EPSG:32608'

import os
import glob

# This code converts velocity fields from some initial CRS into UTM Zone 8N. It
# transforms both the vectors and the underlying pixel coordinates.

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


def projectVelocity(vxFile, vyFile, vvFile, dt):
    '''
    Projects the velocity vectors into the UTM coordinate system.

    Parameters
    ----------
    vxFile, vyFile, vvFile : filenames of the x- and y-components of the 
        velocity in the original reference frame as well as the speed
    dt : time interval used for calculating velocities [a]

    Returns
    -------
    vx_utm, vy_utm, vv_utm : arrays containing the x- and y- components of the 
        velocity and the speed in UTM Zone 8N; the vectors have been projected
        into UTM but the arrays containing the vectors are not projected until
        running writeUTMgeotiff
    '''
    
    # FIX THIS LATER
    # For some strange reason this doesn't seem to affect things too much. Why?
    # dt = 11/365.25 # velocities are in m a^{-1} and timespan between scenes is 11 days
    
    # Import coordinates of pixel centers (X0, Y0), arrays containing velocity
    # components and speed, and the CRS of the original data set.
    X0, Y0, vx, src_crs = importTerraSARX(vxFile)
    _, _, vy, _ = importTerraSARX(vyFile)    
    _, _, vv, _ = importTerraSARX(vvFile)    
    
   
    
    # Create a transformer to project from original coordinates to UTM Zone 8N.
    transformProjection = Transformer.from_crs(src_crs, dst_crs)
    
    # Calculate new coordinates in the original CRS.
    X1 = X0 + vx*dt
    Y1 = Y0 + vy*dt
    
    # Initial and final position in UTM Zone 8N.
    X0_utm, Y0_utm = transformProjection.transform(X0,Y0)
    X1_utm, Y1_utm = transformProjection.transform(X1,Y1)
    
    # Set no data values to NaN because they get turned into 'inf' during transformation.
    X0_utm[np.isnan(vx)], Y0_utm[np.isnan(vx)] = np.nan, np.nan 
    X1_utm[np.isnan(vx)], Y1_utm[np.isnan(vx)] = np.nan, np.nan 
    
    # Compute velocity components and speed.
    vx_utm = (X1_utm-X0_utm)/dt
    vy_utm = (Y1_utm-Y0_utm)/dt
    vv_utm = np.sqrt(vx_utm**2 + vy_utm**2)
    
    return(vx_utm, vy_utm, vv_utm, vx, vy, vv)


def writeUTMgeotiff(data, originalFile):
    '''
    Writes new geotiff file in UTM coordinates, with velocity vectors 
    rotated into the UTM reference frame.

    Parameters
    ----------
    data : array of velocity component or magnitude projected in UTM
        coordinates (the array is not projected, but the velocity is)
    originalFile : Name of the original file; needed in order to extract 
        reference frame information.

    Returns
    -------
    None.
    '''
    
    nx = data.shape[0]
    ny = data.shape[1]
    
    ds = gdal.Open(originalFile)
    dst_ds = gdal.GetDriverByName('GTiff').Create('tmp.tif', ny, nx, 1, gdal.GDT_Float32) # create the 1-band raster file
    dst_ds.SetGeoTransform(ds.GetGeoTransform())    # specify coords and pixel size
    dst_ds.SetProjection(ds.GetSpatialRef().ExportToWkt()) # specify reference frame
    dst_ds.GetRasterBand(1).WriteArray(data)   # write data to the raster
    dst_ds.FlushCache()                     # write to disk
    dst_ds = None
    
    # now warp and crop
    minX=524000
    maxX=604000
    minY=6450000
    maxY=6515000
    
    ## Should be able to run gdal from within python, but for some reason the 
    ## following produces a blank geotiff. Instead use subprocess for now.
    # options = gdal.WarpOptions(srcSRS=ds.GetSpatialRef(), dstSRS='epsg:32608', format='GTiff', xRes=100, yRes=-100, outputBounds=[minX, minY, maxX, maxY])
    # gdal.Warp('test2.tif', 'test.tif', options=options)
    
    fileNameUTM = originalFile[:-4]+'_utm.tif'
    subprocess.run(['gdalwarp', '-t_srs', dst_crs, '-te', str(minX), str(minY), str(maxX), str(maxY), '-tr','100','-100', 'tmp.tif', fileNameUTM])
    os.remove('tmp.tif')

#%%

filenameVX=sorted(glob.glob('/home/jason/Dropbox/Taku/Data/TSX/Release2/Release2/Alaska-Taku/*/*vx_v04.0.tif'))
filenameVY=sorted(glob.glob('/home/jason/Dropbox/Taku/Data/TSX/Release2/Release2/Alaska-Taku/*/*vy_v04.0.tif'))
filenameVV=sorted(glob.glob('/home/jason/Dropbox/Taku/Data/TSX/Release2/Release2/Alaska-Taku/*/*vv_v04.0.tif'))

for j in np.arange(0, len(filenameVX)):
    
    vx_utm, vy_utm, vv_utm, vx, vy, vv = projectVelocity(filenameVX[j], filenameVY[j], filenameVV[j], 11/365.25)
    
    writeUTMgeotiff(vx_utm, filenameVX[j])
    writeUTMgeotiff(vy_utm, filenameVY[j])
    writeUTMgeotiff(vv_utm, filenameVV[j])
