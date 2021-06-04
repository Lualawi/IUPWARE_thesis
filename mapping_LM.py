#! usr/bin/env python 

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import floor


#import data
realisation = 'r1i1p1f2'
model = 'ukesm'
lumip = xr.open_dataset('../CMIP/tasLut_'+ model + '_' + realisation + '_obs.nc')
ds2 = xr.open_dataset('../CMIP/tas_'+ model + '_' + realisation + '_obs.nc')
obs = xr.open_dataset('../CMIP/Tair_merge.nc')

luh = xr.open_dataset('../LAF/LAF_remap_new.nc')
      
#luh needs same lon and lat values as obs and lumip because we used coarsen
#lumip = lumip.sel(lat = lumip.lat[0:-1])
#obs = obs.sel(lat = obs.lat[0:-1])
luh['lat'] = lumip.lat
luh['lon'] = lumip.lon


tasmapped = 0*lumip.tasLut.isel(landuse = 1).copy()#initialize mapped data with zeroes
j = 0
k = 0
for i in range(len(lumip.time.values)):
    k = floor(i/12)
    lumip1 = lumip.copy()
    
    result = tasmapped.isel(time =i)
    luh1 = luh.lut.isel(time = k)

    prim = lumip1.tasLut.isel(landuse =0, time = i)
    crop = lumip1.tasLut.isel(landuse =2, time = i)
    urban = lumip1.tasLut.isel(landuse =3, time = i)
    
    result= result.where(luh1 > 1, prim)
    result= result.where((luh1 > 2) | (luh1 < 2), crop)
    result= result.where((luh1 > 3) | (luh1 < 3), urban)
    
    tasmapped[i,:,:] = result


tasmapped1 = 0*tasmapped.copy() #initialize to trim edges where there is no obs
tasunmapped = ds2.tas

for i in range(len(tasmapped.time.values)):
    tasmapped1[i,:,:] = tasmapped.isel(time =i).where(obs.Tair[i,:,:]>0) #trim edges
    tasunmapped[i,:,:] = tasunmapped.isel(time =i).where(obs.Tair[i,:,:]>0)

#save mapped data into output
mappedtemp = tasmapped1.to_dataset()
unmappedtemp = tasunmapped.to_dataset()

mappedtemp.to_netcdf('../output/remaplaftas_global_'+ model + '.nc')