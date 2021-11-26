#! usr/bin/env python 

import xarray as xr
import pandas as pd
from math import *

ar6 = pd.read_csv('../AR6.csv', index_col = 'names')
ar6_names = list(ar6.index.values)
ar6_names.remove('GIC')

model = 'ukesm'
realisation = 'r1i1p1f2'

rmse_dif = xr.open_dataset('../output/rmse_dif_global_'+model+'.nc')
rmse_difcor = xr.open_dataset('../output/rmse_difcor_global_'+model+'.nc')
rmse_dif_seas = xr.open_dataset('../output/rmse_dif_global_seas_'+model+'.nc')
rmse_difcor_seas = xr.open_dataset('../output/rmse_difcor_global_seas_'+model+'.nc')

kge_dif = xr.open_dataset('../output/kge_dif_global_'+model + '.nc')
kge_difcor = xr.open_dataset('../output/kge_difcor_global_'+model + '.nc')
kge_dif_seas = xr.open_dataset('../output/kge_dif_global_seas_'+model+'.nc')
kge_difcor_seas = xr.open_dataset('../output/kge_dif_global_seas_'+model+'.nc')


ds2 = xr.open_dataset('../CMIP/Tair_merge.nc')


rmses_dif_clim = []
rmses_difcor_clim = []
rmses_dif_jja = []
rmses_dif_son = []
rmses_dif_djf = []
rmses_dif_mam = []

rmses_difcor_jja = []
rmses_difcor_son = []
rmses_difcor_djf = []
rmses_difcor_mam = []

kges_dif_clim = []
kges_difcor_clim = []
kges_dif_jja = []
kges_dif_son = []
kges_dif_djf = []
kges_dif_mam = []

kges_difcor_jja = []
kges_difcor_son = []
kges_difcor_djf = []
kges_difcor_mam = []

for i in ar6_names:
    
    current = i
    lonmin = ar6['Lonmin'][current]
    lonmax = ar6['Lonmax'][current]
    latmin = ar6['Latmin'][current]
    latmax = ar6['Latmax'][current]

    rmsedif = rmse_dif.rmse.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    rmsedifcor = rmse_difcor.rmse.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    rmsedifseas = rmse_dif_seas.rmse.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    rmsedifcorseas = rmse_difcor_seas.rmse.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))

    kgedif = kge_dif.kge.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    kgedifcor = kge_difcor.kge.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    kgedifseas = kge_dif_seas.kge.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))
    kgedifcorseas = kge_difcor_seas.kge.sel(lon = slice(lonmin,lonmax), lat = slice(latmax, latmin))

    rmsedif1 = rmsedif.mean()
    rmsedifcor1 = rmsedifcor.mean()
    
    rmse_djf = rmsedifseas.isel(season = 0).mean()
    rmse_jja = rmsedifseas.isel(season = 1).mean()
    rmse_mam = rmsedifseas.isel(season = 2).mean()
    rmse_son = rmsedifseas.isel(season = 3).mean()
    
    rmse_djfcor = rmsedifcorseas.isel(season = 0).mean()
    rmse_jjacor = rmsedifcorseas.isel(season = 1).mean()
    rmse_mamcor = rmsedifcorseas.isel(season = 2).mean()
    rmse_soncor = rmsedifcorseas.isel(season = 3).mean()
    
    rmses_dif_clim.append(float(rmsedif1))
    rmses_difcor_clim.append(float(rmsedifcor1))
    
    rmses_dif_jja.append(float(rmse_jja))
    rmses_dif_son.append(float(rmse_son))
    rmses_dif_mam.append(float(rmse_mam))
    rmses_dif_djf.append(float(rmse_djf))
    
    rmses_difcor_jja.append(float(rmse_jjacor))
    rmses_difcor_son.append(float(rmse_soncor))
    rmses_difcor_mam.append(float(rmse_mamcor))
    rmses_difcor_djf.append(float(rmse_djfcor))

    kgedif1 = kgedif.mean()
    kgedifcor1 = kgedifcor.mean()
    
    kge_djf = kgedifseas.isel(season = 0).mean()
    kge_jja = kgedifseas.isel(season = 1).mean()
    kge_mam = kgedifseas.isel(season = 2).mean()
    kge_son = kgedifseas.isel(season = 3).mean()
    
    kge_djfcor = kgedifcorseas.isel(season = 0).mean()
    kge_jjacor = kgedifcorseas.isel(season = 1).mean()
    kge_mamcor = kgedifcorseas.isel(season = 2).mean()
    kge_soncor = kgedifcorseas.isel(season = 3).mean()

    kges_dif_clim.append(float(kgedif1))
    kges_difcor_clim.append(float(kgedifcor1))
    
    kges_dif_jja.append(float(kge_jja))
    kges_dif_son.append(float(kge_son))
    kges_dif_mam.append(float(kge_mam))
    kges_dif_djf.append(float(kge_djf))
    
    kges_difcor_jja.append(float(kge_jjacor))
    kges_difcor_son.append(float(kge_soncor))
    kges_difcor_mam.append(float(kge_mamcor))
    kges_difcor_djf.append(float(kge_djfcor))
    
    
rmses_dif_clim = pd.Series(rmses_dif_clim) 
rmses_dif_clim.to_csv('../output/rmse_dif_regional_clim_'+model+'.csv',index = False)
rmses_dif_jja = pd.Series(rmses_dif_jja)
rmses_dif_jja.to_csv('../output/rmse_dif_regional_jja_'+model+'.csv',index = False)
rmses_dif_son = pd.Series(rmses_dif_son) 
rmses_dif_son.to_csv('../output/rmse_dif_regional_son_'+model+'.csv',index = False)
rmses_dif_djf = pd.Series(rmses_dif_djf) 
rmses_dif_djf.to_csv('../output/rmse_dif_regional_djf_'+model+'.csv',index = False)
rmses_dif_mam = pd.Series(rmses_dif_mam) 
rmses_dif_mam.to_csv('../output/rmse_dif_regional_mam_'+model+'.csv',index = False)

rmses_difcor_jja = pd.Series(rmses_difcor_jja)
rmses_difcor_jja.to_csv('../output/rmse_difcor_regional_jja_'+model+'.csv',index = False)
rmses_difcor_son = pd.Series(rmses_difcor_son) 
rmses_difcor_son.to_csv('../output/rmse_difcor_regional_son_'+model+'.csv',index = False)
rmses_difcor_djf = pd.Series(rmses_difcor_djf) 
rmses_difcor_djf.to_csv('../output/rmse_difcor_regional_djf_'+model+'.csv',index = False)
rmses_difcor_mam = pd.Series(rmses_difcor_mam) 
rmses_difcor_mam.to_csv('../output/rmse_difcor_regional_mam_'+model+'.csv',index = False)

rmses_difcor_clim = pd.Series(rmses_difcor_clim)
rmses_difcor_clim.to_csv('../output/rmse_difcor_regional_clim_'+model+'.csv',index = False)

kge_dif = pd.Series(kges_dif_clim) 
kge_dif.to_csv('../output/kge_dif_regional_'+model+'.csv',index = False)
kge_difcor = pd.Series(kges_difcor_clim) 
kge_difcor.to_csv('../output/kge_difcor_regional_'+model+'.csv',index = False)
kges_dif_jja = pd.Series(kges_dif_jja)
kges_dif_jja.to_csv('../output/kge_dif_regional_jja_'+model+'.csv',index = False)
kges_dif_son = pd.Series(kges_dif_son) 
kges_dif_son.to_csv('../output/kge_dif_regional_son_'+model+'.csv',index = False)
kges_dif_djf = pd.Series(kges_dif_djf) 
kges_dif_djf.to_csv('../output/kge_dif_regional_djf_'+model+'.csv',index = False)
kges_dif_mam = pd.Series(kges_dif_mam) 
kges_dif_mam.to_csv('../output/kge_dif_regional_mam_'+model+'.csv',index = False)

kges_difcor_jja = pd.Series(kges_dif_jja)
kges_difcor_jja.to_csv('../output/kge_difcor_regional_jja_'+model+'.csv',index = False)
kges_difcor_son = pd.Series(kges_difcor_son) 
kges_difcor_son.to_csv('../output/kge_difcor_regional_son_'+model+'.csv',index = False)
kges_difcor_djf = pd.Series(kges_difcor_djf) 
kges_difcor_djf.to_csv('../output/kge_difcor_regional_djf_'+model+'.csv',index = False)
kges_difcor_mam = pd.Series(kges_difcor_mam) 
kges_difcor_mam.to_csv('../output/kge_difcor_regional_mam_'+model+'.csv',index = False)