    

import os
import numpy as np
import rasterio
import glob
import re
from datetime import datetime
import datetime as dt
from statsmodels.nonparametric.smoothers_lowess import lowess

def stack_geotiffs_and_find_max_index(directory, date_structure, filename_ending,output_directory, year,datatype, mask_switch):
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
    # Create an array that has the number of observations (excluding NaNs)
    # Find the maximum value for each cell across the 3D array
    max_values = np.nanmax(stacked_arrays, axis=0)
    min_values = np.nanmin(stacked_arrays, axis=0)
    velocity_range = max_values - min_values
    num_observations1 = np.sum(~np.isnan(stacked_arrays), axis=0)

    # Identify the date corresponding to the maximum value for each cell
    max_indices = np.argmax(stacked_arrays, axis=0)
    max_dates = np.array(sorted_dates)[max_indices]
    # Convert max_dates to day of year
    max_dates_doy = np.vectorize(lambda date: datetime.strptime(date, '%Y-%m-%d').timetuple().tm_yday)(max_dates)
    
    if mask_switch == 'True':
     # Filter max_dates to be between day 32 and day 274
        mask = (max_dates_doy >= 32) & (max_dates_doy <= 274)
        max_values_filtered = np.where(mask, max_dates_doy, np.nan)
        max_dates_doy = max_values_filtered

    # Save the velocity_range and day of year array as a GeoTIFF file
    output_filepath1 = os.path.join(output_directory, 'day_of_max_'f'{year}'f'_{datatype}_redo.tif')
    output_filepath2 = os.path.join(output_directory, 'velocity_range'f'{year}'f'_{datatype}_redo.tif')
    output_filepath3 = os.path.join(output_directory, 'num_obs'f'{year}'f'_{datatype}_redo.tif')
    meta.update(dtype=rasterio.float32, count=1)
    with rasterio.open(output_filepath1, 'w', **meta) as dst:
        dst.write(max_dates_doy.astype(rasterio.float32), 1)
    with rasterio.open(output_filepath2, 'w', **meta) as dst:
        dst.write(velocity_range.astype(rasterio.float32), 1)
    print(f"Saved max values array to {output_filepath1}")
    with rasterio.open(output_filepath3, 'w', **meta) as dst:
        dst.write(num_observations1.astype(rasterio.float32), 1)
    print(f"Saved max values array to {output_filepath1}")

    return stacked_arrays, sorted_dates, max_values, max_dates, max_dates_doy


if __name__ == '__main__':

    datatype = 'TSX'
    if datatype == 'Sentinel1':
        date_structure = r's1cycle_(\d{2}[A-Za-z]{3}\d{2})_(\d{2}[A-Za-z]{3}\d{2})_'
        directory = '/hdd3/taku/satellite_velocity/Sentinel1/Release1-12day/'
        filename_ending = 'vv_v02.0_utm.tif'
    if datatype == 'TSX':
        date_structure = r'_(\d{2}[A-Za-z]{3}\d{2})_(\d{2}[A-Za-z]{3}\d{2})_'
        directory = '/hdd/taku/TerraSAR-X/velocities/Release4/Alaska-Taku/'
        filename_ending = 'vv_v04.0_utm.tif'

    mask_switch = 'True'
    year = '2022'
    output_directory = ''
    stacked_arrays, sorted_dates, max_values, max_dates, max_dates_doy = stack_geotiffs_and_find_max_index(directory, date_structure, filename_ending,output_directory, year, datatype, mask_switch)

    dates = [dt.datetime.strptime(j,'%Y-%m-%d') for j in sorted_dates]
    ordinal_dates = [x.toordinal() for x in dates]

j = 300
k = 300

lowess(dates, stacked_arrays[:,j,k], frac=0.1, it=3)

