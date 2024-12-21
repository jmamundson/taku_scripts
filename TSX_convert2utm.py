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
from os.path import basename

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

#%%
def rotateVelocity(vxFile, vyFile):
    '''
    Rotates the velocity vector from polar stereographic coordinates to UTM. 
    This requires finding the angle of the position vector in polar 
    stereographic coordinates and rotating the coordinate system so that one
    velocity component points north-south and the other points east west. 
    
    Note that in polar stereographic coordinates, for the orientation used for
    the Taku imagery, negative y points roughly from the pole toward the equator
    and x points roughly from west to east. These directions are accounted for
    when determining the UTM velocity, so that vx_utm points east and vy_utm
    points north.
    
    Parameters
    ----------
    vxFile, vyFile : filenames of the x- and y-components of the velocity in 
            the original reference frame
     
    Returns
    -------
    vx_utm, vy_utm : components of the velocity vector in UTM coordinates, with
            positive values indicative motion to the east and north, 
            respectively

    '''
    
    X0, Y0, vx, src_crs = importTerraSARX(vxFile)
    _, _, vy, _ = importTerraSARX(vyFile)  

    vx_utm = np.zeros(vx.shape)
    vy_utm = np.zeros(vy.shape)

    # x-coordinate points from west to east
    # negative y-coordinate points from north to south (starting at north pole)

    theta = np.arctan(X0/-Y0) # rotation angle
    
    # there might be a more compact way of stepping through every pixel using matrix multiplication        
    for j in np.arange(0, theta.shape[0]):
        for k in np.arange(0, theta.shape[1]):
            rot = np.array([ [np.cos(theta[j,k]), np.sin(theta[j,k])] , [-np.sin(theta[j,k]), np.cos(theta[j,k])] ])        
            vy_utm[j,k], vx_utm[j,k] = np.matmul(rot, np.array([-vy[j,k], vx[j,k]]))

    vy_utm = -vy_utm
    
    
    return(vx_utm, vy_utm)

#%%


# def projectVelocity(vxFile, vyFile, dt):
#     '''
#     Projects the velocity vectors into the UTM coordinate system.

#     Parameters
#     ----------
#     vxFile, vyFile : filenames of the x- and y-components of the 
#         velocity in the original reference frame
#     dt : time interval used for calculating velocities [a]

#     Returns
#     -------
#     vx_utm, vy_utm, vv_utm : arrays containing the x- and y- components of the 
#         velocity and the speed in UTM Zone 8N; the vectors have been projected
#         into UTM but the arrays containing the vectors are not projected until
#         running writeUTMgeotiff
#     '''
    
        
#     # Import coordinates of pixel centers (X0, Y0), arrays containing velocity
#     # components and speed, and the CRS of the original data set.
#     X0, Y0, vx, src_crs = importTerraSARX(vxFile)
#     _, _, vy, _ = importTerraSARX(vyFile)        
   
    
#     # Create a transformer to project from original coordinates to UTM Zone 8N.
#     transformProjection = Transformer.from_crs(src_crs, dst_crs)
    
#     # Calculate new coordinates in the original CRS.
#     X1 = X0 + vx*dt
#     Y1 = Y0 + vy*dt
    
#     # Initial and final position in UTM Zone 8N.
#     X0_utm, Y0_utm = transformProjection.transform(X0,Y0)
#     X1_utm, Y1_utm = transformProjection.transform(X1,Y1)
    
#     # Set no data values to NaN because they get turned into 'inf' during transformation.
#     X0_utm[np.isnan(vx)], Y0_utm[np.isnan(vx)] = np.nan, np.nan 
#     X1_utm[np.isnan(vx)], Y1_utm[np.isnan(vx)] = np.nan, np.nan 
    
#     # Compute velocity components and speed.
#     vx_utm = (X1_utm-X0_utm)/dt
#     vy_utm = (Y1_utm-Y0_utm)/dt
        
#     return(vx_utm, vy_utm)


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
    
    filenameUTM = originalFile[:-4]+'_utm.tif'
    subprocess.run(['gdalwarp', '-t_srs', dst_crs, '-te', str(minX), str(minY), str(maxX), str(maxY), '-tr','100','-100', '-r', 'cubicspline', 'tmp.tif', filenameUTM])
    os.remove('tmp.tif')

    return(filenameUTM)


def writeUTMgeotiff_vv(data, filename):
    '''
    Create new geotiff. Does not gdalwarp, so use this only for creating a 
    vv geotiff, after creating vx and vy geotiffs.

    Parameters
    ----------
    data : array of speed
    filename : name of the original file; needed in order to create new filename

    Returns
    -------
    None.

    '''
    
    filenameUTM = filename[:-4]+'_utm.tif'
    
    driver = gdal.GetDriverByName('GTiff')
    rows, cols = data.shape
    dataset = driver.Create(filenameUTM, cols, rows, 1, gdal.GDT_Float32)
    
    # now warp and crop
    minX=524000
    maxX=604000
    minY=6450000
    maxY=6515000
    xres = 100
    yres = 100


    geotransform = (minX, xres, 0, maxY, 0, -yres)  # Example geotransform
    projection = osr.SpatialReference()
    projection.ImportFromEPSG(32608)  # Example projection (WGS84, UTM Zone 8)
    
    # Set the geotransform and projection
    dataset.SetGeoTransform(geotransform)
    dataset.SetProjection(projection.ExportToWkt())

    # Write the array data to the GeoTIFF file
    band = dataset.GetRasterBand(1)
    band.WriteArray(data)
    band.FlushCache()

    # Close the dataset
    dataset = None
    
    

#%%

datatype = 'Sentinel1'

if datatype == 'TSX':
    base_dir = '/hdd/taku/TerraSAR-X/velocities/Release4/Alaska-Taku/'
    filenameVX=sorted(glob.glob(base_dir + '*/*vx_v04.0.tif'))
    filenameVY=sorted(glob.glob(base_dir + '*/*vy_v04.0.tif'))
    filenameVV=sorted(glob.glob(base_dir + '*/*vv_v04.0.tif'))
elif datatype == 'Sentinel1':
    base_dir = '/hdd/taku/Sentinel1/Release1-12day/'
    filenameVX=sorted(glob.glob(base_dir + '*/*vx_v02.0.tif'))
    filenameVY=sorted(glob.glob(base_dir + '*/*vy_v02.0.tif'))
    filenameVV=sorted(glob.glob(base_dir + '*/*vv_v02.0.tif'))

# dt = 11/365.25

for j in np.arange(0, len(filenameVX)):
    
    print('Processing: ' + basename(filenameVX[j]))
    # velocity vectors are in UTM, but pixel coordinates are not
    # vx_utm, vy_utm = projectVelocity(filenameVX[j], filenameVY[j], dt) 
    vx_utm, vy_utm = rotateVelocity(filenameVX[j], filenameVY[j])
    
    
    # project pixel coordinates into UTM and save geotiffs
    filenameVX_new = writeUTMgeotiff(vx_utm, filenameVX[j])
    filenameVY_new = writeUTMgeotiff(vy_utm, filenameVY[j])
    
    # load new geotiffs and compute the speed; now vectors and pixel coordinates are already in UTM
    X, Y, vx, src_crs = importTerraSARX(filenameVX_new)
    X, Y, vy, src_crs = importTerraSARX(filenameVY_new)
    
    vv = np.sqrt(vx**2+vy**2)
    
    writeUTMgeotiff_vv(vv, filenameVV[j])

