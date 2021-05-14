#! usr/bin/env python 

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import floor

#import data
realisation = 'r1i1p1f1'
model = 'cesm2'
lumip = xr.open_dataset('../CMIP/tasLut_'+ model + '_' + realisation + '_obs.nc')
ds2 = xr.open_dataset('../CMIP/tas_'+ model + '_' + realisation + '_obs.nc')
obs = xr.open_dataset('../CMIP/Tair_merge.nc')

luh = xr.open_dataset('../LAF/LAF_obs_new.nc')

      
#luh needs same lon and lat values as obs and lumip because we used coarsen
#lumip = lumip.sel(lat = lumip.lat[0:-1])
#obs = obs.sel(lat = obs.lat[0:-1])
luh['lat'] = lumip.lat
luh['lon'] = lumip.lon

#check if summation of landuse types amounts to 1
check = luh.prim + luh.urban + luh.crop #+ luh.pastr

tasmapped = 0*lumip.tasLut.isel(landuse= 0).copy()#initialize mapped data with zeroes
j = 0
k = 0
for i in range(len(lumip.time.values)):
    k = floor(i/12)
    lumip1 = lumip.copy()
    prim = lumip1.tasLut.isel(landuse=0,time =i)*luh.prim.isel(time =k)
    ds = xr.Dataset({'prim': prim} )
    #ds['pastr'] = lumip1.tasLut.isel(landuse=1, time = i)*luh.pastr.isel(time =k)
    ds['crop'] = lumip1.tasLut.isel(landuse=2, time = i)*luh.crop.isel(time =k)
    
    ds['urban'] = lumip1.tasLut.isel(landuse=3, time=i)*luh.urban.isel(time =k)
    
    vars_to_sum = ['prim', 'crop', 'urban'] #'pastr',
    result = ds[vars_to_sum].to_array().sum("variable") #skip na is already true
    
    #here you need to remove points where urban area exists but lumip is empty
    #before adding the resulting tas values otherwise
    # you add partial data even coz urban area takes portion without data
    zer = 0*result.copy()
    urban = luh.urban.isel(time =k)
    tasurban = lumip1.tasLut.isel(landuse=3, time=i)
    zer = zer.where((urban>0) & (np.isnan(tasurban)), 1)
    result = result * zer
    result = result.where(result>0) #remove values multiplied by zero
    
    #same for crop and prim
    zer1 = 0*result.copy()
    prim1 = luh.prim.isel(time =k)
    tasprim = lumip1.tasLut.isel(landuse=0, time=i)
    zer1 = zer1.where((prim1>0) & (np.isnan(tasprim)), 1)
    result = result * zer1
    result = result.where(result>0)
    
    zer2 = 0*result.copy()
    crops = luh.crop.isel(time =k)
    tascrop = lumip1.tasLut.isel(landuse=2, time=i)
    zer2 = zer2.where((crops>0) & (np.isnan(tascrop)), 1)
    result = result * zer2
    result = result.where(result>0)
    
    tasmapped[i,:,:] = result


tasmapped1 = 0*tasmapped.copy() #initialize to trim edges where there is no obs
tasunmapped = ds2.tas

for i in range(len(tasmapped.time.values)):
    tasmapped1[i,:,:] = tasmapped.isel(time =i).where(obs.Tair[i,:,:]>0) #trim edges
    tasunmapped[i,:,:] = tasunmapped.isel(time =i).where(obs.Tair[i,:,:]>0)

#save mapped data into output
mappedtemp = tasmapped1.to_dataset()
unmappedtemp = tasunmapped.to_dataset()

mappedtemp.to_netcdf('../output/mappedtas_global_'+ model + '.nc')
unmappedtemp.to_netcdf('../output/unmappedtas_global_'+ model + '.nc')
    

