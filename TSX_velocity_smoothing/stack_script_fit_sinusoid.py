import os
import numpy as np
import rasterio
import glob
import re
from datetime import datetime
import datetime as dt
from statsmodels.nonparametric.smoothers_lowess import lowess
from osgeo import gdal
import subprocess
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import fmin, minimize, LinearConstraint
from matplotlib import pyplot as plt


#%%
f = 1/365.25 # frequency [1/yr], expressed in days
omega = 2*np.pi*f # angular frequency [1/yr], expressed in days
daysYear = 365.25

def findPhi(X, t, v):
    A, phi = X
    
    res = np.sqrt(np.sum((-A*np.cos(omega*t + phi) - v)**2))

    return(res)


#%%
def stack_geotiffs_and_fit_sinusoid(directory, date_structure, filename_ending,output_directory, datatype, mask_switch, gaussian_switch):
    # Get a list of all subdirectories
    subdirectories = next(os.walk(directory))[1]

    # Sort the list just in case the directories are not in the order you expect
    subdirectories.sort()

    arrays_dict = {} # load arrays into this dict, then stack them later

    # Iterate over the subdirectories
    for subdir in subdirectories:
        # Skip directories named 'npz' and 'plots'
        if subdir in ['npz', 'plots', 'empty_files']:
            print(f"Skipping directory: {subdir}")
            continue

        # Find all files that end with the specified filename ending
        geotiff_files = glob.glob(f'{directory}{subdir}/*{filename_ending}')
        print(geotiff_files)
        
        # Print the filenames for debugging
        if not geotiff_files:
            print(f"No GeoTIFF files found in {subdir} with ending {filename_ending}")
            continue
        
        # Process only the first GeoTIFF file found
        filepath = geotiff_files[0]
        
        
        filename = os.path.basename(filepath)
        match = re.search(date_structure, filename)
        if match:
            date_str = match.group(1)
            date_obj = datetime.strptime(date_str, '%d%b%y')
            formatted_date = date_obj.strftime('%Y-%m-%d')
            print(f"Processing file: {filename}, Date: {formatted_date}")

        with rasterio.open(filepath) as src:
            band1 = src.read(1)
            print('Band1 has shape', band1.shape)
            arrays_dict[formatted_date] = band1
            meta = src.meta.copy()  # Copy the metadata for later use

    # Stack the arrays into a 3D array
    sorted_dates = sorted(arrays_dict.keys())
    stacked_arrays = np.stack([arrays_dict[date] for date in sorted_dates], axis=0)
    
        
    # express dates as ordinal dates
    dates = np.array([dt.datetime.strptime(j,'%Y-%m-%d') for j in sorted_dates])
    ordinal_dates = [x.toordinal() for x in dates]
    
        
    # spatially smooth the velocity fields
    if gaussian_switch == 'True':
        for j in np.arange(0,len(dates)):
            stacked_arrays[j,:,:] = gaussian_filter(stacked_arrays[j,:,:],0.1,order=0)
            
       
            
    phi = np.zeros(band1.shape)
    A = np.zeros(band1.shape)
    u_mean = np.zeros(band1.shape)
    
    totalPixels = A.shape[0]*A.shape[1]
    
    if datatype == 'Sentinel1':
        t_mask = (dates>dt.datetime(2016,12,31)) & (dates<dt.datetime(2022,1,1))
        jmin = 120
        jmax = 210
        kmin = 110
        kmax = 180
        
    elif datatype == 'TSX':
        t_mask = (dates>dt.datetime(2020,12,31)) & (dates<dt.datetime(2025,1,1))
        jmin = 240
        jmax = 410
        kmin = 220
        kmax = 360
    
    # this little bit fits a sinusoid through every grid point
    for j in np.arange(0,stacked_arrays.shape[1]):
        for k in np.arange(0,stacked_arrays.shape[2]):
            if (j>jmin) and (j<jmax): # don't run on full data set
                if (k>kmin) and (k<kmax):
                    t = np.array(ordinal_dates)
                    v = stacked_arrays[:,j,k]
                    
                    # apply mask to only analyze data during specific time period
                    v = v[t_mask]
                    t = t[t_mask]
                    
                    # remove nan values
                    t = t[~np.isnan(v)]
                    v = v[~np.isnan(v)]
                    
                    # subtract mean value during time period
                    u_mean[j,k] = np.mean(v)
                    v = v-np.mean(v)
                    
                                                            
                    # minimization, with constraints on A and phi                    
                    cons = LinearConstraint([[1,0],[0,1]], [0,-np.pi], [np.inf, np.pi])
                    result = minimize(findPhi, (100,0), args=(t,v), constraints=[cons])
                    
                    A[j,k], phi[j,k] = result.x
                    print("{:.2f}".format( (j-jmin)*(k-kmin) / ((jmax-jmin)*(kmax-kmin)) ) )
                    
    # convert phi to day of year
    doy = daysYear/2 * (1 - phi/np.pi)
    
    
    # Save the velocity_range and day of year array as a GeoTIFF file
    
    output_filepath1 = os.path.join(output_directory, datatype + '_doy.tif')
    output_filepath2 = os.path.join(output_directory, datatype + '_A.tif')
    output_filepath3 = os.path.join(output_directory, datatype + '_uMean.tif')
        
    
    takuShp = 'taku_outline_large_poly.shp'
    
    meta.update(dtype=rasterio.float32, count=1)
    
    if os.path.exists('*.tif'):
        os.remove('*.tif')
    
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(doy.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath1])
    os.remove('tmp.tif')
    print(f"Saved doy array to {output_filepath1}")
        
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(A.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath2])
    os.remove('tmp.tif')
    print(f"Saved amplitude values array to {output_filepath2}")
    
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(u_mean.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath3])
    os.remove('tmp.tif')
    print(f"Saved u_mean values array to {output_filepath3}")
    


    return(stacked_arrays, dates, phi, A, u_mean)


    

#%%
if __name__ == '__main__':

    datatype = 'Sentinel1'
    # datatype = 'TSX'
    if datatype == 'Sentinel1':
        date_structure = r's1cycle_(\d{2}[A-Za-z]{3}\d{2})_(\d{2}[A-Za-z]{3}\d{2})_'
        directory = '/hdd/taku/Sentinel1/Release1-12day/'
        filename_ending = 'vv_v02.0_utm.tif'
    if datatype == 'TSX':
        date_structure = r'_(\d{2}[A-Za-z]{3}\d{2})_(\d{2}[A-Za-z]{3}\d{2})_'
        directory = '/hdd/taku/TerraSAR-X/velocities/Release4/Alaska-Taku/'
        filename_ending = 'vv_v04.0_utm.tif'

    gaussian_switch = 'False'
    mask_switch = 'False'
    output_directory = ''
    
    
    stacked_arrays, dates, phi, A, u_mean = stack_geotiffs_and_fit_sinusoid(directory, date_structure, filename_ending, output_directory, datatype, mask_switch, gaussian_switch)

#%%
with rasterio.open(datatype + '_doy.tif') as src:
    doy = src.read(1)
    
    doy[doy==0] = np.nan

with rasterio.open(datatype + '_A.tif') as src:
    A = src.read(1)
    A[A==0] = np.nan
    
with rasterio.open(datatype + '_uMean.tif') as src:
    u_mean = src.read(1)
    u_mean[A==0] = np.nan
    
#doy[A<1] = np.nan

if datatype == 'Sentinel1':
    ymin = 100
    ymax = 220
    xmin = 110
    xmax = 190
    
elif datatype == 'TSX':
    ymin = 240
    ymax = 410
    xmin = 220
    xmax = 360

doy[A<5] = np.nan
A[A<5] = np.nan


plt.figure(figsize=(10,5))
plt.subplot(131)
plt.imshow(doy,vmin=50,vmax=250)
plt.xlim([xmin, xmax])
plt.ylim([ymax, ymin])
cbar = plt.colorbar(location='bottom')
cbar.set_label('Day of peak speed')

plt.subplot(132)
plt.imshow(A)
plt.xlim([xmin, xmax])
plt.ylim([ymax, ymin])
cbar = plt.colorbar(location='bottom')
cbar.set_label('Amplitude of oscillations [m/yr]')

plt.subplot(133)
plt.imshow(A/u_mean)
plt.xlim([xmin, xmax])
plt.ylim([ymax, ymin])
cbar = plt.colorbar(location='bottom')
cbar.set_label('Amplitude / mean speed')

plt.savefig(datatype + '_phase_fitting.png', format='png', dpi=150)

#%%
plt.figure()
j = 380 # north-south pixel
k = 310 # east-west pixel
ordinal_dates = [x.toordinal() for x in dates]
plt.plot(dates, stacked_arrays[:,j,k], dates, -A[j,k]*np.cos(omega*np.array(ordinal_dates) + phi[j,k]) + u_mean[j,k])