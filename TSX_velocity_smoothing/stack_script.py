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

def stack_geotiffs_and_find_max_index(directory, date_structure, filename_ending,output_directory, year,datatype, mask_switch, gaussian_switch):
    # Get a list of all subdirectories
    subdirectories = next(os.walk(directory))[1]

    # Sort the list just in case the directories are not in the order you expect
    subdirectories.sort()

    # Filter subdirectories to include only those from the year 2022
    subdirectories = [subdir for subdir in subdirectories if year in subdir]
    
    #edit if you only want to process a few subdirectories
    last_three_subdirectories = subdirectories[:]
    arrays_dict = {}

    # Iterate over the subdirectories
    for subdir in last_three_subdirectories:
        # Skip directories named 'npz' and 'plots'
        if subdir in ['npz', 'plots', 'empty_files']:
            print(f"Skipping directory: {subdir}")
            continue

        # Find all files that end with the specified filename ending
        geotiff_files = glob.glob(f'{directory}{subdir}/*{filename_ending}')
        
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
    dates = [dt.datetime.strptime(j,'%Y-%m-%d') for j in sorted_dates]
    ordinal_dates = [x.toordinal() for x in dates]
    
    if gaussian_switch == 'True':
        for j in np.arange(0,len(dates)):
            stacked_arrays[j,:,:] = gaussian_filter(stacked_arrays[j,:,:],0.1,order=0)
            
            
            
    # filter the stacked_arrays using a lowess filter
    data = np.zeros(stacked_arrays.shape)
    for j in np.arange(0,stacked_arrays.shape[1]):
        for k in np.arange(0,stacked_arrays.shape[2]):
            data[:,j,k] = lowess(stacked_arrays[:,j,k], ordinal_dates,  frac=1/4, it=1, return_sorted=False)
    
    stacked_arrays = data            
            
    # find when max velocity occurs
    max_indices = np.argmax(stacked_arrays, axis=0)
    
    #max_dates = dates[max_indices]
    max_dates_doy = np.array(ordinal_dates)[max_indices]
    max_dates_doy -= dt.datetime(int(year)-1,12,31).toordinal() # express max_dates in terms of day of year
   
    max_values = np.nanmax(stacked_arrays, axis=0)
    min_values = np.nanmin(stacked_arrays, axis=0)

    velocity_range = max_values - min_values
    num_observations1 = np.sum(~np.isnan(stacked_arrays), axis=0)

    # # Identify the date corresponding to the maximum value for each cell
    # max_indices = np.argmax(stacked_arrays, axis=0)
    # max_dates = np.array(sorted_dates)[max_indices]
    # # Convert max_dates to day of year
    # max_dates_doy = np.vectorize(lambda date: datetime.strptime(date, '%Y-%m-%d').timetuple().tm_yday)(max_dates)
    

    
    if mask_switch == 'True':
     # Filter max_dates to be between day 32 and day 274
        mask = (max_dates_doy >= 60) & (max_dates_doy <=305)
        max_dates_doy = np.where(mask, max_dates_doy, np.nan)
        max_values = np.where(mask, max_values, np.nan)
        min_values = np.where(mask, min_values, np.nan)
        
    # Save the velocity_range and day of year array as a GeoTIFF file
    output_filepath1 = os.path.join(output_directory, 'day_of_max_'f'{year}'f'_{datatype}_redo.tif')
    output_filepath2 = os.path.join(output_directory, 'max_values_'f'{year}'f'_{datatype}_redo.tif')
    output_filepath3 = os.path.join(output_directory, 'min_values_'f'{year}'f'_{datatype}_redo.tif')
    output_filepath4 = os.path.join(output_directory, 'num_obs_'f'{year}'f'_{datatype}_redo.tif')
    
    
    takuShp = 'taku_outline_large_poly.shp'
    
    meta.update(dtype=rasterio.float32, count=1)
    
    if os.path.exists('*' + year + '*.tif'):
        os.remove('*.tif')
    
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(max_dates_doy.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath1])
    os.remove('tmp.tif')
    print(f"Saved max dates array to {output_filepath1}")
        
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(max_values.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath2])
    os.remove('tmp.tif')
    print(f"Saved max values array to {output_filepath1}")
    
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(min_values.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath3])
    os.remove('tmp.tif')
    print(f"Saved min values array to {output_filepath1}")
    
    with rasterio.open('tmp.tif', 'w', **meta) as dst:
        dst.write(num_observations1.astype(rasterio.float32), 1)
    subprocess.run(['gdalwarp', '-cutline', takuShp, '-dstnodata', '0', 'tmp.tif', output_filepath4])
    os.remove('tmp.tif')
    print(f"Saved num observations array to {output_filepath1}")


    return stacked_arrays, sorted_dates, max_values, max_dates_doy, min_values


    


if __name__ == '__main__':

    datatype = 'Sentinel1'
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
    year = '2018'
    output_directory = ''
    stacked_arrays, sorted_dates, max_values, max_dates_doy, min_values = stack_geotiffs_and_find_max_index(directory, date_structure, filename_ending, output_directory, year, datatype, mask_switch, gaussian_switch)

#%%
from matplotlib import pyplot as plt
with rasterio.open('day_of_max_2018_Sentinel1_redo.tif') as src:
    day_max = src.read(1)
with rasterio.open('num_obs_2018_Sentinel1_redo.tif') as src:
    number_obs = src.read(1)
    
#day_max[number_obs<10] = np.nan    
plt.imshow(day_max,vmin=90,vmax=365-90)
plt.colorbar()