#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt

model = 'ukesm'
realisation = 'r1i1p1f2'

obs = xr.open_dataset('../CMIP/Tair_merge.nc')
umpd = xr.open_dataset('../output/unmappedtas_global_' + model + '.nc')
mpd = xr.open_dataset('../output/mappedtas_global_' + model + '.nc')
umpdcor = xr.open_dataset('../output/tasunmappedcor_global_' + model + '.nc')
mpdcor = xr.open_dataset('../output/tasmappedcor_global_' + model + '.nc')

umpd['time'] = obs.time
mpd['time'] = obs.time
umpdcor['time'] = obs.time
mpdcor['time'] = obs.time

obsclim = obs.Tair.mean(dim = 'time')
umpdclim = umpd.tas.mean(dim = 'time')
mpdclim = mpd.tasLut.mean(dim = 'time')
umpdcorclim = umpdcor.tas.mean(dim = 'time')
mpdcorclim = mpdcor.tasLut.mean(dim = 'time')

obsseas = obs.Tair.groupby('time.season').mean('time')
umpdseas = umpd.tas.groupby('time.season').mean('time')
mpdseas = mpd.tasLut.groupby('time.season').mean('time')
umpdcorseas = umpdcor.tas.groupby('time.season').mean('time')
mpdcorseas = mpdcor.tasLut.groupby('time.season').mean('time')

rmse_umpd = ((((umpd.tas - obs.Tair)**2).mean('time'))**0.5).rename('rmse')
rmse_mpd = ((((mpd.tasLut - obs.Tair)**2).mean('time'))**0.5).rename('rmse')
rmse_umpdcor = ((((umpdcor.tas - obs.Tair)**2).mean('time'))**0.5).rename('rmse')
rmse_mpdcor = ((((mpdcor.tasLut - obs.Tair)**2).mean('time'))**0.5).rename('rmse')

rmse_umpd.to_netcdf('../output/rmse_umpd_global_'+model+'.nc')
rmse_mpd.to_netcdf('../output/rmse_mpd_global_'+model+'.nc')
rmse_umpdcor.to_netcdf('../output/rmse_umpdcor_global_'+model+'.nc')
rmse_mpdcor.to_netcdf('../output/rmse_mpdcor_global_'+model+'.nc')

rmse_dif = (rmse_mpd - rmse_umpd).rename('rmse')
rmse_difcor = (rmse_mpdcor - rmse_umpdcor).rename('rmse')
rmse_dif.to_netcdf('../output/rmse_dif_global_'+model+'.nc')
rmse_difcor.to_netcdf('../output/rmse_difcor_global_'+model+'.nc')


rmse_umpd_seas = ((((umpd.tas - obs.Tair)**2).groupby('time.season').mean('time'))**0.5).rename('rmse')
rmse_mpd_seas = ((((mpd.tasLut - obs.Tair)**2).groupby('time.season').mean('time'))**0.5).rename('rmse')
rmse_umpdcor_seas = ((((umpdcor.tas - obs.Tair)**2).groupby('time.season').mean('time'))**0.5).rename('rmse')
rmse_mpdcor_seas = ((((mpdcor.tasLut - obs.Tair)**2).groupby('time.season').mean('time'))**0.5).rename('rmse')

rmse_dif = (rmse_mpd_seas - rmse_umpd_seas).rename('rmse')
rmse_difcor = (rmse_mpdcor_seas - rmse_umpdcor_seas).rename('rmse')

rmse_umpd_seas.to_netcdf('../output/rmse_umpd_global_seas_'+model + '.nc')
rmse_mpd_seas.to_netcdf('../output/rmse_mpd_global_seas_'+model + '.nc')
rmse_umpdcor_seas.to_netcdf('../output/rmse_umpdcor_global_seas_'+model + '.nc')
rmse_mpdcor_seas.to_netcdf('../output/rmse_mpdcor_global_seas_'+model + '.nc')

rmse_dif.to_netcdf('../output/rmse_dif_global_seas_'+model+'.nc')
rmse_difcor.to_netcdf('../output/rmse_difcor_global_seas_'+model+'.nc')
    
B_umpd = umpdclim/obsclim
B_mpd = mpdclim/obsclim
B_umpdcor = umpdcorclim/obsclim
B_mpdcor = mpdcorclim/obsclim

umpd['time'] = obs.time
mpd['time'] = obs.time
umpdcor['time'] = obs.time
mpdcor['time'] = obs.time

G_umpd = (umpd.tas.std(dim = 'time')/umpdclim)/(obs.Tair.std(dim = 'time') / obsclim)
G_mpd = (mpd.tasLut.std(dim = 'time')/mpdclim)/(obs.Tair.std(dim = 'time') / obsclim)
G_umpdcor = (umpdcor.tas.std(dim = 'time')/umpdcorclim)/(obs.Tair.std(dim = 'time') / obsclim)
G_mpdcor = (mpdcor.tasLut.std(dim = 'time')/mpdcorclim)/(obs.Tair.std(dim = 'time') / obsclim)


P_umpd = xr.corr(umpd.tas , obs.Tair, dim= 'time')
P_mpd = xr.corr(mpd.tasLut , obs.Tair, dim= 'time')
P_umpdcor = xr.corr(umpdcor.tas , obs.Tair, dim= 'time')
P_mpdcor = xr.corr(mpdcor.tasLut , obs.Tair, dim= 'time')

kge_umpd = (1-((P_umpd-1)**2 + (B_umpd - 1)**2 + (G_umpd-1)**2)**(0.5)).rename('kge')
kge_mpd = (1-((P_mpd-1)**2 + (B_mpd - 1)**2 + (G_mpd-1)**2)**(0.5)).rename('kge')
kge_umpdcor = (1-((P_umpdcor-1)**2 + (B_umpdcor - 1)**2 + (G_umpdcor-1)**2)**(0.5)).rename('kge')
kge_mpdcor = (1-((P_mpdcor-1)**2 + (B_mpdcor - 1)**2 + (G_mpdcor-1)**2)**(0.5)).rename('kge')

kge_umpd.to_netcdf('../output/kge_umpd_global_'+model + '.nc')
kge_mpd.to_netcdf('../output/kge_mpd_global_'+model + '.nc')
kge_umpdcor.to_netcdf('../output/kge_umpdcor_global_'+model + '.nc')
kge_mpdcor.to_netcdf('../output/kge_mpdcor_global_'+model + '.nc')

kge_dif = (kge_mpd - kge_umpd).rename('kge')
kge_difcor = (kge_mpdcor - kge_umpdcor).rename('kge')

kge_dif.to_netcdf('../output/kge_dif_global_'+model + '.nc')
kge_difcor.to_netcdf('../output/kge_difcor_global_'+model + '.nc')

kge_umpd_seas = 0*rmse_umpd_seas.copy().rename('kge')
kge_mpd_seas = 0*rmse_mpd_seas.copy().rename('kge')
kge_umpdcor_seas = 0*rmse_umpdcor_seas.copy().rename('kge')
kge_mpdcor_seas = 0*rmse_mpdcor_seas.copy().rename('kge')
kge_dif_seas = 0*rmse_umpd_seas.copy().rename('kge')
kge_difcor_seas = 0*rmse_mpdcor_seas.copy().rename('kge')

for i in range (4):
    seas = ['DJF','JJA','MAM', 'SON']
    umpd_seas = umpd.sel(time=umpd.time.dt.season==seas[i])
    mpd_seas = mpd.sel(time=mpd.time.dt.season==seas[i])
    umpdcor_seas = umpdcor.sel(time=umpdcor.time.dt.season==seas[i])
    mpdcor_seas = mpdcor.sel(time=mpdcor.time.dt.season==seas[i])
    obs_seas = obs.sel(time=obs.time.dt.season==seas[i])
    
    obsclim_seas = obs_seas.Tair.mean(dim = 'time')
    umpdclim_seas = umpd_seas.tas.mean(dim = 'time')
    mpdclim_seas = mpd_seas.tasLut.mean(dim = 'time')
    umpdcorclim_seas = umpdcor_seas.tas.mean(dim = 'time')
    mpdcorclim_seas = mpdcor_seas.tasLut.mean(dim = 'time')
    
    G_umpd = (umpd_seas.tas.std(dim = 'time')/umpdclim_seas)/(obs_seas.Tair.std(dim = 'time') / obsclim_seas)
    G_mpd = (mpd_seas.tasLut.std(dim = 'time')/mpdclim_seas)/(obs_seas.Tair.std(dim = 'time') / obsclim_seas)
    G_umpdcor = (umpdcor_seas.tas.std(dim = 'time')/umpdcorclim_seas)/(obs_seas.Tair.std(dim = 'time') / obsclim_seas)
    G_mpdcor = (mpdcor_seas.tasLut.std(dim = 'time')/mpdcorclim_seas)/(obs_seas.Tair.std(dim = 'time') / obsclim_seas)
    
    
    P_umpd = xr.corr(umpd_seas.tas , obs_seas.Tair, dim= 'time')
    P_mpd = xr.corr(mpd_seas.tasLut , obs_seas.Tair, dim= 'time')
    P_umpdcor = xr.corr(umpdcor_seas.tas , obs_seas.Tair, dim= 'time')
    P_mpdcor = xr.corr(mpdcor_seas.tasLut , obs_seas.Tair, dim= 'time')
    
    kge_umpd = 1-((P_umpd-1)**2 + (B_umpd - 1)**2 + (G_umpd-1)**2)**(0.5)
    kge_mpd = 1-((P_mpd-1)**2 + (B_mpd - 1)**2 + (G_mpd-1)**2)**(0.5)
    kge_umpdcor = 1-((P_umpdcor-1)**2 + (B_umpdcor - 1)**2 + (G_umpdcor-1)**2)**(0.5)
    kge_mpdcor = 1-((P_mpdcor-1)**2 + (B_mpdcor - 1)**2 + (G_mpdcor-1)**2)**(0.5)

    kge_umpd_seas[i,:,:] = kge_umpd
    kge_mpd_seas[i,:,:] = kge_mpd
    kge_umpdcor_seas[i,:,:] = kge_umpdcor
    kge_mpdcor_seas[i,:,:] = kge_mpdcor
    
    kge_dif = kge_mpd - kge_umpd
    kge_difcor = kge_mpdcor - kge_umpdcor
    
    kge_dif_seas[i,:,:] = kge_dif
    kge_difcor_seas[i,:,:] = kge_difcor

kge_umpd_seas.to_netcdf('../output/kge_umpd_global_seas_'+model + '.nc')
kge_mpd_seas.to_netcdf('../output/kge_mpd_global_seas_'+model + '.nc')
kge_umpdcor_seas.to_netcdf('../output/kge_umpdcor_global_seas_'+model + '.nc')
kge_mpdcor_seas.to_netcdf('../output/kge_mpdcor_global_seas_'+model + '.nc')

kge_dif_seas.to_netcdf('../output/kge_dif_global_seas_'+model + '.nc')
kge_difcor_seas.to_netcdf('../output/kge_difcor_global_seas_'+model + '.nc')
