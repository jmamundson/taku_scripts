#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 13:45:35 2025

@author: jason
"""

import geoutils as gu
import numpy as np
import matplotlib.pyplot as plt
import xdem
import glob
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Coregisters all DEMs in ./02_takuGrid/ to the 20250822 DEM using a vertical shift.")
    
        
    if len(sys.argv) == 2: # one command-line argument has been passed to the script, in addition to the script's own name
        parser.print_help()
        sys.exit(0) 
        
    base = '/hdd2/taku/uav/surveys'
    # Open a reference and to-be-aligned DEM
    ref_dem = xdem.DEM(base + '/02_takuGrid/taku_DEM_20250822.tif')
    ref_dem.save(base + '/02_coregister/' + ref_dem.name.split('/')[-1])

    dems = sorted(glob.glob(base + '/02_takuGrid/taku_DEM_20*.tif'))
    dems = np.append(dems, ['/hdd2/taku/DEMs_orthomosaics/20150802_ChrisLarsen/DEMs/taku_DEM_20150802.tif'])
    mask = base + '/02_coregister/landingStrip.shp'

    for j in np.arange(0,len(dems)):
    	
        if j in [2,7,9]:
            continue
        
        tba_dem = xdem.DEM(dems[j])
    
        if ref_dem.name==tba_dem.name: # don't coregister reference image to itself
    	    print('Skip reference image.')
    	    continue
        
        bedrockShp = gu.Vector(mask)
        bedrockMask = bedrockShp.create_mask(ref_dem) # glacier mask
    
        # mycoreg = xdem.coreg.NuthKaab() # + xdem.coreg.Deramp(poly_order=1)
        # mycoreg = xdem.coreg.Deramp(poly_order=1) #
        mycoreg = xdem.coreg.VerticalShift()
    
        mycoreg.fit(ref_dem, tba_dem, inlier_mask=bedrockMask)
        dem_aligned = mycoreg.apply(tba_dem)
        
        # elevation difference
        # dh = ref_dem - dem_aligned
    
    
        dem_aligned.save(base + '/02_coregister/' + dems[j].split('/')[-1]
                         )
    
   
