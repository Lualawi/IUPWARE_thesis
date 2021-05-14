#! usr/bin/env python 

import xarray as xr
import numpy as np
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from scipy import stats


model = 'cesm2'
realisation = 'r1i1p1f1'
tas = xr.open_dataset('../CMIP/tas_'+model+'_180.nc')
tas = tas.sel(lat = slice(-56,84))

obs = xr.open_dataset('../CMIP/Tair_merge.nc')
obs = obs.sel(lat = slice(84,-56))

#import CESM data for grid reference
elev_sub = xr.open_dataset('../elev/elev_obs.nc')
elev_grid = xr.open_dataset('../elev/elev_grid_' + model + '.nc')

remap = 'mapped'
tasmapped = xr.open_dataset('../output/'+remap+'tas_global_' + model + '.nc')
tasmapped = tasmapped.sel(lat = slice(84,-56))
tasunmapped = xr.open_dataset('../output/unmappedtas_global_'+ model +'.nc')
tasunmapped = tasunmapped.sel(lat = slice(84,-56))

#determine the number of lon and lat at GCM level inside study area
reglon = len(tas.lon.values)
reglat = len(tas.lat.values)

#initialize lapse rate
lapse = 0*obs.Tair
reg = 0*obs.Tair

for j in range(0,reglon-1,1):
    for k in range (0,reglat-1,1):
        #slice the data into 1 grid size of GCM model

        lonslice = slice(tas.lon.values[j],tas.lon.values[j+1])

        latslice = slice(tas.lat.values[k+1],tas.lat.values[k])
            
        obslice = obs.sel(lon = lonslice, lat = latslice)
        elevslice = elev_sub.sel(lon = lonslice, lat =latslice)
        lapslice = lapse.sel(lon = lonslice, lat = latslice)
        #regslice = reg.sel(lon = lonslice, lat = latslice)
        
        #stack dataset
        stacked_obs = obslice.Tair.stack(z = ("lat","lon"))
        stacked_elev = elevslice.Band1.stack(z = ("lat","lon")).values
        stacked_lapse = lapslice.stack(z = ("lat","lon"))
        #stacked_reg  = regslice.stack(z = ("lat","lon"))
        
        # apply regression between observation and elevation for each timestep and determine slope
        for l in range(144):
            stacked_obs1 = stacked_obs.isel(time =l).values
            mask = ~np.isnan(stacked_elev) & ~np.isnan(stacked_obs1)
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(stacked_elev[mask],stacked_obs1[mask])
                #rsq = r_value**2
            except:
                slope = np.nan
                #r_value = np.nan
            stacked_lapse.isel(time = l).values.fill(slope)
            #stacked_reg.isel(time = l).values.fill(rsq)
            
        #unstacked lapse data for lon lat that is filled with the slope from regression
        unstacked_lapse = stacked_lapse.unstack("z")
        #unstacked_reg = stacked_reg.unstack("z")
        lapse.loc[dict(lon = lonslice, lat = latslice)] = unstacked_lapse
        #reg.loc[dict(lon = lonslice, lat = latslice)] = unstacked_reg


#now define the correction by multiplying lapse rate with elevation diference
elev_difer = elev_sub.Band1 - elev_grid.Band1
correction = elev_difer * lapse


tascor = tasmapped.tasLut.copy(deep =True)
tascor1 = tasunmapped.tas.copy(deep = True)

for m in range(len(tasmapped.time.values)):
        tascor[m,:,:] = tasmapped.tasLut.isel(time = m) + correction.isel(time = m)
        tascor1[m,:,:] = tasunmapped.tas.isel(time = m) + correction.isel(time =m)
        tascor[m,:,:] = tascor.isel(time = m).where(obs.Tair.isel(time = m) >0)
        tascor1[m,:,:] = tascor1.isel(time = m).where(obs.Tair.isel(time = m) >0)

#    rsq_list[i] = reg.mean()
tascor.to_netcdf('../output/tas'+remap+'cor_global_'+ model + '.nc')
tascor1.to_netcdf('../output/tasunmappedcor_global_' + model + '.nc')


