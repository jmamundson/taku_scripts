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


def projectVelocity(vxFile, vyFile, dt):
    '''
    Projects the velocity vectors into the UTM coordinate system.

    Parameters
    ----------
    vxFile, vyFile, vvFile : filenames of the x- and y-components of the 
        velocity in the original reference frame
    dt : time interval used for calculating velocities [a]

    Returns
    -------
    vx_utm, vy_utm, vv_utm : arrays containing the x- and y- components of the 
        velocity and the speed in UTM Zone 8N; the vectors have been projected
        into UTM but the arrays containing the vectors are not projected until
        running writeUTMgeotiff
    '''
    
        
    # Import coordinates of pixel centers (X0, Y0), arrays containing velocity
    # components and speed, and the CRS of the original data set.
    X0, Y0, vx, src_crs = importTerraSARX(vxFile)
    _, _, vy, _ = importTerraSARX(vyFile)        
   
    
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
        
    return(vx_utm, vy_utm)


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
base_dir = '/hdd/taku/TerraSAR-X/velocities'
filenameVX=sorted(glob.glob(base_dir + '/Release4/Alaska-Taku/*/*vx_v04.0.tif'))
filenameVY=sorted(glob.glob(base_dir + '/Release4/Alaska-Taku/*/*vy_v04.0.tif'))
filenameVV=sorted(glob.glob(base_dir + '/Release4/Alaska-Taku/*/*vv_v04.0.tif'))


dt = 11/365.25

for j in np.arange(0, len(filenameVX)):
    
    # velocity vectors are in UTM, but pixel coordinates are not
    vx_utm, vy_utm = projectVelocity(filenameVX[j], filenameVY[j], dt) 
    
    # project pixel coordinates into UTM and save geotiffs
    filenameVX_new = writeUTMgeotiff(vx_utm, filenameVX[j])
    filenameVY_new = writeUTMgeotiff(vy_utm, filenameVY[j])
    
    # load new geotiffs and compute the speed; now vectors and pixel coordinates are already in UTM
    X, Y, vx, src_crs = importTerraSARX(filenameVX_new)
    X, Y, vy, src_crs = importTerraSARX(filenameVY_new)
    
    vv = np.sqrt(vx**2+vy**2)
    
    writeUTMgeotiff_vv(vv, filenameVV[j])



#%%
# from shapely.geometry import LineString, Point
# import shapefile
# import rasterstats
# from rasterstats import zonal_stats, point_query


# profile_shp = '/home/jason/Downloads/cross_transect_test.shp'
# profile = shapefile.Reader(profile_shp).shapes()

# profile = LineString(profile[0].points)
# X, Y = np.array(profile.coords.xy)
# X, Y = np.linspace(X[0], X[-1], 100, endpoint=True), np.linspace(Y[0], Y[-1], 100, endpoint=True) 

# profile = LineString(list(zip(X,Y)))

# # profiles_shp = '../sediment_profiles/profiles.shp'
# # profiles = shapefile.Reader(profiles_shp).shapes()
# # profile = LineString(profiles[3].points)
# # X, Y = np.array(profile.coords.xy)


# filenameVX = sorted(glob.glob('/home/jason/Desktop/Sentinel1/Release1-12day/Vel-2023-07-30.2023-08-10/*vx_v02.0_utm.tif'))
# filenameVY = sorted(glob.glob('/home/jason/Desktop/Sentinel1/Release1-12day/Vel-2023-07-30.2023-08-10/*vy_v02.0_utm.tif'))
# filenameVV = sorted(glob.glob('/home/jason/Desktop/Sentinel1/Release1-12day/Vel-2023-07-30.2023-08-10/*vv_v02.0_utm.tif'))

# # ds = gdal.Open(filenameVX[0])
# # vx = ds.ReadAsArray()

# # ds = gdal.Open(filenameVY[0])
# # vy = ds.ReadAsArray()

# # ds = gdal.Open(filenameVV[0])
# # vv = ds.ReadAsArray()

# # plt.imshow(np.sqrt(vx**2+vy**2)-vv)
# # plt.colorbar()


# dist = np.sqrt((X-X[0])**2+(Y-Y[0])**2)
# vx = np.array(rasterstats.point_query(profile, filenameVX[0])[0])
# vy = np.array(rasterstats.point_query(profile, filenameVY[0])[0])
# vv = np.array(rasterstats.point_query(profile, filenameVV[0])[0])

# plt.subplot(211)
# plt.plot(dist, np.sqrt(vx**2+vy**2), dist, vv)
# plt.subplot(212)
# plt.plot(dist, np.sqrt(vx**2+vy**2)-vv)
