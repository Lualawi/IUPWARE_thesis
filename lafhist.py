#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import seaborn as sns
from scipy import stats

sns.set_style("white")

model1 = 'ukesm'
model2 = 'cesm2'
realisation = 'r1i1p1f2'

rmse_dif1 = xr.open_dataset('../output/rmse_dif_laf_'+model1+'.nc')
rmse_difcor1 = xr.open_dataset('../output/rmse_difcor_laf_'+model1+'.nc')
kge_delta1 = xr.open_dataset('../output/kge_dif_laf_'+model1 + '.nc')
kge_deltacor1 = xr.open_dataset('../output/kge_difcor_laf_'+model1 + '.nc')

rmse_umpdcor1= xr.open_dataset('../output/rmse_umpdcor_laf_'+model1+'.nc')
rmse_mpdcor1= xr.open_dataset('../output/rmse_mpdcor_global_'+model1+'.nc')
kge_umpdcor1= xr.open_dataset('../output/kge_umpd_laf_'+model1 + '.nc')
kge_mpdcor1= xr.open_dataset('../output/kge_mpd_global_'+model1 + '.nc')

rmse_dif2 = xr.open_dataset('../output/rmse_dif_laf_'+model2+'.nc')
rmse_difcor2 = xr.open_dataset('../output/rmse_difcor_laf_'+model2+'.nc')
kge_delta2 = xr.open_dataset('../output/kge_dif_laf_'+model2 + '.nc')
kge_deltacor2 = xr.open_dataset('../output/kge_difcor_laf_'+model2 + '.nc')

rmse_umpdcor2= xr.open_dataset('../output/rmse_umpdcor_laf_'+model2+'.nc')
rmse_mpdcor2= xr.open_dataset('../output/rmse_mpdcor_global_'+model2+'.nc')
kge_umpdcor2= xr.open_dataset('../output/kge_umpd_laf_'+model2 + '.nc')
kge_mpdcor2= xr.open_dataset('../output/kge_mpd_global_'+model2 + '.nc')



obs = xr.open_dataset('../CMIP/Tair_merge.nc')
obsseas = obs.Tair.groupby('time.season').mean('time')

cmap_whole = plt.cm.get_cmap('RdBu_r')
cmap55 = cmap_whole(0.01)
cmap50 = cmap_whole(0.05)   #blue
cmap45 = cmap_whole(0.1)
cmap40 = cmap_whole(0.15)
cmap35 = cmap_whole(0.2)
cmap30 = cmap_whole(0.25)
cmap25 = cmap_whole(0.3)
cmap20 = cmap_whole(0.325)
cmap10 = cmap_whole(0.4)
cmap5 = cmap_whole(0.45)
cmap1 = cmap_whole(0.49)
cmap_1 = cmap_whole(0.51)
cmap0 = cmap_whole(0.5)
cmap_5 = cmap_whole(0.55)
cmap_10 = cmap_whole(0.6)
cmap_20 = cmap_whole(0.625)
cmap_25 = cmap_whole(0.7)
cmap_30 = cmap_whole(0.75)
cmap_35 = cmap_whole(0.8)
cmap_40 = cmap_whole(0.85)
cmap_45 = cmap_whole(0.9)
cmap_50 = cmap_whole(0.95)  #red
cmap_55 = cmap_whole(0.99)


colors = [cmap45,cmap30,cmap10, cmap1, cmap_1,cmap_10,cmap_30,cmap_45]
colors1 = [cmap_45,cmap_30,cmap_10, cmap_1, cmap1,cmap10,cmap30,cmap45]
    
    # declare list of colors for discrete colormap of colorbar
cmap = mpl.colors.ListedColormap(colors,N=len(colors))
cmap.set_over(cmap_55)
cmap.set_under(cmap55)
cmap.set_bad(color='0.8')

cmap2 = mpl.colors.ListedColormap(colors1,N=len(colors1))
cmap2.set_over(cmap55)
cmap2.set_under(cmap_55)
cmap.set_bad(color='0.8')

# colorbar args
values1 = [-1.5,-1.0,-0.5,-0.25,0,0.25,0.5,1.0,1.5]
values2 = [-0.15,-0.10,-0.05,-0.025,0,0.025,0.05,0.10,0.15]
tick_locs1 = [-1.5,-1.0,-0.5,0,0.5,1.0,1.5]
tick_labels1 = ['-1.5','-1.0','-0.5','0','0.5','1.0','1.5']
tick_locs2 = ['-0.15','-0.10','-0.05','0','0.05','0.10','0.15']
norm1 = mpl.colors.BoundaryNorm(values1,cmap.N)
norm2 = mpl.colors.BoundaryNorm(values2,cmap2.N)

rmse_dif1 = rmse_dif1.reindex(lat=list(reversed(rmse_dif1.lat)))
rmse_difcor1 = rmse_difcor1.reindex(lat=list(reversed(rmse_difcor1.lat)))

kge_delta1 = kge_delta1.reindex(lat=list(reversed(kge_delta1.lat)))
kge_deltacor1 = kge_deltacor1.reindex(lat=list(reversed(kge_deltacor1.lat)))

rmse_umpdcor1 = rmse_umpdcor1.rmse.stack(z = ("lat","lon"))
rmse_mpdcor1 = rmse_mpdcor1.rmse.stack(z = ("lat","lon"))
kge_umpdcor1 = kge_umpdcor1.kge.stack(z = ("lat","lon"))
kge_mpdcor1 = kge_mpdcor1.kge.stack(z = ("lat","lon"))
    
lats = rmse_dif1.lat.values
lons = rmse_dif1.lon.values
lons, lats = np.meshgrid(lons,lats)

lats1 = rmse_difcor1.lat.values
lons1 = rmse_difcor1.lon.values
lons1, lats1 = np.meshgrid(lons1,lats1)

rmse_dif2 = rmse_dif2.reindex(lat=list(reversed(rmse_dif2.lat)))
rmse_difcor2 = rmse_difcor2.reindex(lat=list(reversed(rmse_difcor2.lat)))

kge_delta2 = kge_delta2.reindex(lat=list(reversed(kge_delta2.lat)))
kge_deltacor2 = kge_deltacor2.reindex(lat=list(reversed(kge_deltacor2.lat)))

rmse_umpdcor2 = rmse_umpdcor2.rmse.stack(z = ("lat","lon"))
rmse_mpdcor2 = rmse_mpdcor2.rmse.stack(z = ("lat","lon"))
kge_umpdcor2 = kge_umpdcor2.kge.stack(z = ("lat","lon"))
kge_mpdcor2 = kge_mpdcor2.kge.stack(z = ("lat","lon"))
    
lats = rmse_dif2.lat.values
lons = rmse_dif2.lon.values
lons, lats = np.meshgrid(lons,lats)

lats2 = rmse_difcor2.lat.values
lons2 = rmse_difcor2.lon.values
lons2, lats2 = np.meshgrid(lons2,lats2)



num_bins = np.arange(0,20,0.5)
num_bins2 = np.arange(-0.41,1,0.05)

a = rmse_umpdcor1.values[~np.isnan(rmse_umpdcor1.values)]
b = rmse_mpdcor1.values[~np.isnan(rmse_mpdcor1.values)]
c = kge_umpdcor1.values[~np.isnan(kge_umpdcor1.values)]
d = kge_mpdcor1.values[~np.isnan(kge_mpdcor1.values)]

e = rmse_umpdcor2.values[~np.isnan(rmse_umpdcor2.values)]
f = rmse_mpdcor2.values[~np.isnan(rmse_mpdcor2.values)]
g = kge_umpdcor2.values[~np.isnan(kge_umpdcor2.values)]
h = kge_mpdcor2.values[~np.isnan(kge_mpdcor2.values)]


kde1 = stats.gaussian_kde(a)
kde2 = stats.gaussian_kde(b)
x1 = np.linspace(0, 20, 1000)

kde3 = stats.gaussian_kde(e)
kde4 = stats.gaussian_kde(f)
x2 = np.linspace(0, 20, 1000)

kde5 = stats.gaussian_kde(c)
kde6 = stats.gaussian_kde(d)
x3 = np.linspace(-0.41, 1, 1000)

kde7 = stats.gaussian_kde(g)
kde8 = stats.gaussian_kde(h)
x4 = np.linspace(-0.41, 1, 1000)
fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_subplot(221)
ax1.spines['bottom'].set_color('0.8')
ax1.spines['top'].set_color('0.8')
ax1.spines['left'].set_color('0.8')
ax1.spines['right'].set_color('0.8')
ax1.plot(x1, kde1(x1), color = 'blue', label = 'LM')
ax1.plot(x1, kde2(x1), color = 'red', label = 'WM')
ax1.legend()
ax1.set_ylabel( 'RMSE')
ax1.set_title('UKESM1',fontweight = 'bold')
ax1.set_title('a.',fontweight = 'bold', loc='left')

ax2 = fig.add_subplot(222)
ax2.spines['bottom'].set_color('0.8')
ax2.spines['top'].set_color('0.8')
ax2.spines['left'].set_color('0.8')
ax2.spines['right'].set_color('0.8')
ax2.plot(x2, kde3(x2), color = 'blue', label = 'LM')
ax2.plot(x2, kde4(x2), color = 'red', label = 'WM')
ax2.legend()
ax2.set_title('CESM2',fontweight = 'bold')
ax2.set_title('b.',fontweight = 'bold', loc='left')

ax3 = fig.add_subplot(223)
ax3.spines['bottom'].set_color('0.8')
ax3.spines['top'].set_color('0.8')
ax3.spines['left'].set_color('0.8')
ax3.spines['right'].set_color('0.8')
ax3.plot(x3, kde5(x3), color = 'blue', label = 'LM')
ax3.plot(x3, kde6(x3), color = 'red', label = 'WM')

ax3.legend()
ax3.set_ylabel( 'KGE')
ax3.set_title('c.',fontweight = 'bold', loc='left')

ax4 = fig.add_subplot(224)
ax4.spines['bottom'].set_color('0.8')
ax4.spines['top'].set_color('0.8')
ax4.spines['left'].set_color('0.8')
ax4.spines['right'].set_color('0.8')
ax4.plot(x4, kde7(x4), color = 'blue', label = 'LM')
ax4.plot(x4, kde8(x4), color = 'red', label = 'WM')

ax4.legend()
ax4.set_title('d.',fontweight = 'bold', loc='left')

plt.savefig('../plot_jun/metrics_laf_kde_both.png')






