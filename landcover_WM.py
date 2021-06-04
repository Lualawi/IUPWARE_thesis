#! /usr/bin/env python

import xarray as xr
import pandas as pd
import numpy as np

years = np.arange(2001,2015)
#AR6 regions
ar6 = pd.read_csv('../AR6.csv', index_col = 'names')
ar6_names = list(ar6.index.values)


for i in years:
  name = '../ESA/ESACCI-LC-L4-LCCS-Map-300m-P1Y-'+str(i)+'-v2.0.7b.nc'
  
  #import landuse
  landuse = xr.open_dataset(name)
  
  lu_classes = landuse.lccs_class
  #    classes = pd.DataFrame(data = lu_classes.flag_meanings.split(), index = lu_classes.flag_values)
      
      
  #reclassify land cover values to lumip tiles (4 tiles)
  land1 = lu_classes.where(lu_classes > 0)
  land2 = land1.where((land1 != 10) & (land1!=20) & (land1 !=30), 3)
  land3 = land2.where((land2 > 180) | (land2 <10), 1)
  land4 = land3.where(land3 != 190, 4)
  land5 = land4.where((land4 <200) | (land4>209),1)
  land6 = land5.where((land5 <205) & (land5 >0))
  
  
  crop = land6.where(land6 == 3,0)
  crop = crop.where(land6 != 3,1)
  crop = crop.coarsen(lat=18, boundary="trim").sum()
  crop = crop.coarsen(lon=18, boundary="trim").sum()
  
  prim = land6.where(land6 == 1,0)
  prim = prim.coarsen(lat=18, boundary="trim").sum()
  prim = prim.coarsen(lon=18, boundary="trim").sum()
  
 # pasture = land7.where(land7 == 2,0)
 # pasture = pasture.where(land7 != 2,1)
 # pasture = pasture.coarsen(lat=18, boundary="trim").sum()
 # pasture = pasture.coarsen(lon=18, boundary="trim").sum()
  
  urban = land6.where(land6 == 4,0)
  urban = urban.where(land6 != 4,1)
  urban = urban.coarsen(lat=18, boundary="trim").sum()
  urban = urban.coarsen(lon=18, boundary="trim").sum()
  
  landsum = crop + prim + urban
  
  laf = crop/landsum
  laf = laf.rename('crop')
  laf = laf.to_dataset()
  laf['prim'] = prim/landsum
 # laf['pastr'] = pasture/landsum
  laf['urban'] = urban/landsum
  laf = laf.assign_coords(time=np.datetime64(str(i)))
  laf = laf.expand_dims('time')
  
  laf.to_netcdf('../LAF/'+'LAF_'+str(i)+'_coarsened.nc')