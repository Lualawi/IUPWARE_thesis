#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl

model = 'ukesm'
realisation = 'r1i1p1f2'

rmse_dif = xr.open_dataset('../output/rmse_dif_global_'+model+'.nc')
rmse_difcor = xr.open_dataset('../output/rmse_difcor_global_'+model+'.nc')
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


colors = [cmap45,cmap40,cmap30,cmap20,cmap10, cmap1, cmap_1,cmap_10,cmap_20,cmap_30,cmap_40,cmap_45]
colors1 = [cmap_45,cmap_40,cmap_30,cmap_20,cmap_10, cmap_1, cmap1,cmap10,cmap20,cmap30,cmap40,cmap45]
    
    # declare list of colors for discrete colormap of colorbar
cmap = mpl.colors.ListedColormap(colors,N=len(colors))
cmap.set_over(cmap_55)
cmap.set_under(cmap55)
cmap.set_bad(color='0.8')

cmap2 = mpl.colors.ListedColormap(colors1,N=len(colors1))
cmap2.set_over(cmap55)
cmap2.set_under(cmap_55)
cmap2.set_bad(color='0.8')

# colorbar args
values1 = [-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5]
values2 = [-0.15,-0.125,-0.10,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.10,0.125,0.15]
tick_locs1 = [-1.5,-1.0,-0.5,0,0.5,1.0,1.5]
tick_labels1 = ['-1.5','-1.25','-1.0','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1.0','1.25','1.5']
tick_labels2 = ['-0.15','-0.125','-0.10','-0.075','-0.05','-0.025','0','0.025','0.05','0.075','0.10','0.125','0.15']
norm1 = mpl.colors.BoundaryNorm(values1,cmap.N)
norm2 = mpl.colors.BoundaryNorm(values2,cmap2.N)

rmse_dif = rmse_dif.reindex(lat=list(reversed(rmse_dif.lat)))
rmse_difcor = rmse_difcor.reindex(lat=list(reversed(rmse_difcor.lat)))

kge_delta = xr.open_dataset('../output/kge_dif_global_'+model + '.nc')
kge_deltacor = xr.open_dataset('../output/kge_difcor_global_'+model + '.nc')

kge_delta = kge_delta.reindex(lat=list(reversed(kge_delta.lat)))
kge_deltacor = kge_deltacor.reindex(lat=list(reversed(kge_deltacor.lat)))
    

lats = rmse_dif.lat.values
lons = rmse_dif.lon.values
lons, lats = np.meshgrid(lons,lats)

lats1 = rmse_difcor.lat.values
lons1 = rmse_difcor.lon.values
lons1, lats1 = np.meshgrid(lons1,lats1)

fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_axes([0.05,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,rmse_dif.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_title('Root Mean Squared Error (not corrected)')
ax1.set_title('a.',fontweight = 'bold', loc='left')

cbax = fig.add_axes([0.05, 0.2, 0.45, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{mapped}$ - $RMSE_{unmapped}$', fontsize=9)
cb.ax.set_xticklabels(tick_labels1)
cb.ax.tick_params(labelsize=7)
#plot arrows
bluelabel = 'Less error/Higher Skill'
redlabel = 'More error/Lower Skill'

cb.ax.text(0.7, -3, redlabel, size=8, ha='center', va='center',transform  = cbax.transAxes)
cb.ax.text(0.3, -3, bluelabel, size=8, ha='center', va='center', transform  = cbax.transAxes)

cb.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=redlabel,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)
cb.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=bluelabel,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)

ax2 = fig.add_axes([0.53,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax2
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,kge_delta.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax2.set_title('Kling-Gupta Efficiency (not corrected)')
ax2.set_title('b.',fontweight = 'bold', loc='left')

cbax2 = fig.add_axes([0.53, 0.2, 0.45, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{mapped}$ - $KGE_{unmapped}$', fontsize=9)
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)
#plot arrows
bluelabel1 = 'Higher efficiency/Higher Skill'
redlabel1 = 'Lower efficiency/Lower Skill'

cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax2.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax2.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)


plt.savefig('../plot1/global/metrics_difmap_nocor_'+model+'.png')

fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_axes([0.05,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,rmse_difcor.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_title('Root Mean Squared Error')
ax1.set_title('a.',fontweight = 'bold', loc='left')

cbax = fig.add_axes([0.05, 0.2, 0.45, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{mapped}$ - $RMSE_{unmapped}$', fontsize=9)
cb.ax.set_xticklabels(tick_labels1)
cb.ax.tick_params(labelsize=7)
#plot arrows
bluelabel = 'Less error/Higher Skill'
redlabel = 'More error/Lower Skill'

cb.ax.text(0.7, -3, redlabel, size=8, ha='center', va='center',transform  = cbax.transAxes)
cb.ax.text(0.3, -3, bluelabel, size=8, ha='center', va='center', transform  = cbax.transAxes)

cb.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=redlabel,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)
cb.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=bluelabel,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)


ax2 = fig.add_axes([0.53,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax2
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,kge_deltacor.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax2.set_title('Kling-Gupta Efficiency')
ax2.set_title('b.',fontweight = 'bold', loc='left')

cbax2 = fig.add_axes([0.53, 0.2, 0.45, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{mapped}$ - $KGE_{unmapped}$', fontsize=9)
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)
#plot arrows
bluelabel1 = 'Higher efficiency/Higher Skill'
redlabel1 = 'Lower efficiency/Lower Skill'

cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax2.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax2.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)

plt.savefig('../plot1/global/metrics_difmap_cor_'+model+'.png')


lon_seas = []
rmse_difcor_seas = xr.open_dataset('../output/rmse_difcor_global_seas_'+model+'.nc')
rmse_difcor_seas = rmse_difcor_seas.reindex(lat=list(reversed(rmse_difcor_seas.lat)))

kge_difcor_seas = xr.open_dataset('../output/kge_difcor_global_seas_'+model+'.nc')
kge_difcor_seas = kge_difcor_seas.reindex(lat=list(reversed(kge_difcor_seas.lat)))


fig = plt.figure(figsize =(10, 7))
ax = fig.add_axes([0.05,0.05,0.9,0.9])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

for i in range(4):
    alp = ['a.', 'b.', 'c.', 'd.']
    pos = [[0.1,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.1,0.13,0.4,0.4],[0.55,0.13,0.4,0.4]]
    
    rmse_dif_si = rmse_difcor_seas.isel(season = i)
    
    a = rmse_dif_si.rmse.mean('lon')
    lon_seas.append(a)
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution=None)
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons1,lats1,rmse_dif_si.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
    ax1.set_title( str(rmse_difcor_seas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax = fig.add_axes([0.25, 0.1, 0.5, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE  = $RMSE_{mapped}$ - $RMSE_{unmapped}$')
cb.ax.set_xticklabels(tick_labels1)
cb.ax.tick_params(labelsize=7)
#plot arrows

cb.ax.text(0.7, -3, redlabel, size=8, ha='center', va='center',transform  = cbax.transAxes)
cb.ax.text(0.3, -3, bluelabel, size=8, ha='center', va='center', transform  = cbax.transAxes)

cb.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=redlabel,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)
cb.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=bluelabel,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax.transAxes)

plt.savefig('../plot1/global/rmse_difmap_cor_season_'+model+'.png') 


fig = plt.figure(figsize =(10, 7))
ax = fig.add_axes([0.05,0.05,0.9,0.9])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

lon_kg_seas = []

for i in range(4):
    alp = ['a.', 'b.', 'c.', 'd.']
    pos = [[0.1,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.1,0.13,0.4,0.4],[0.55,0.13,0.4,0.4]]
    
    kge_dif_si = kge_difcor_seas.isel(season = i)
    
    a = kge_dif_si.kge.mean('lon')
    lon_kg_seas.append(a)
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution=None)
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons1,lats1,kge_dif_si.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
    ax1.set_title( str(kge_difcor_seas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax2 = fig.add_axes([0.25, 0.1, 0.5, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE  = $KGE_{mapped}$ - $KGE_{unmapped}$')
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)
#plot arrows
bluelabel1 = 'Higher efficiency/Higher Skill'
redlabel1 = 'Lower efficiency/Lower Skill'

cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax2.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax2.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)

plt.savefig('../plot1/global/kge_difmap_cor_season_'+model+'.png') 

rmse_umpd = xr.open_dataset('../output/rmse_umpd_global_'+model+'.nc')
rmse_mpd= xr.open_dataset('../output/rmse_mpd_global_'+model+'.nc')
rmse_umpdcor= xr.open_dataset('../output/rmse_umpdcor_global_'+model+'.nc')
rmse_mpdcor= xr.open_dataset('../output/rmse_mpdcor_global_'+model+'.nc')

kge_umpd= xr.open_dataset('../output/kge_umpd_global_'+model + '.nc')
kge_mpd= xr.open_dataset('../output/kge_mpd_global_'+model + '.nc')
kge_umpdcor= xr.open_dataset('../output/kge_umpdcor_global_'+model + '.nc')
kge_mpdcor= xr.open_dataset('../output/kge_mpdcor_global_'+model + '.nc')

rmse1 = rmse_umpd.rmse.stack(z = ("lat","lon"))
rmse2 = rmse_mpd.rmse.stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor.rmse.stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor.rmse.stack(z = ("lat","lon"))
rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

data_rmse = [rmse1, rmse2, rmse3, rmse4]

#plot box and longitude for KGE
kge1 = kge_umpd.kge.stack(z = ("lat","lon"))
kge2 = kge_mpd.kge.stack(z = ("lat","lon"))
kge3 = kge_umpdcor.kge.stack(z = ("lat","lon"))
kge4 = kge_mpdcor.kge.stack(z = ("lat","lon"))
kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

data_kge = [kge1, kge2, kge3, kge4]

rmse_lon_umpd = rmse_umpd.rmse.mean(dim = 'lon')
rmse_lon_mpd = rmse_mpd.rmse.mean(dim = 'lon')
rmse_lon_dif = (rmse_lon_mpd - rmse_lon_umpd).rename('rmse')

rmse_lon_umpdcor = rmse_umpdcor.rmse.mean(dim = 'lon')
rmse_lon_mpdcor = rmse_mpdcor.rmse.mean(dim = 'lon')
rmse_lon_difcor = (rmse_lon_mpdcor - rmse_lon_umpdcor).rename('rmse')

kge_lon_umpd = kge_umpd.kge.mean(dim = 'lon')
kge_lon_mpd = kge_mpd.kge.mean(dim = 'lon')
kge_lon_dif = (kge_lon_mpd - kge_lon_umpd).rename('kge')

kge_lon_umpdcor = kge_umpdcor.kge.mean(dim = 'lon')
kge_lon_mpdcor = kge_mpdcor.kge.mean(dim = 'lon')
kge_lon_difcor = (kge_lon_mpdcor - kge_lon_umpdcor).rename('kge')

#zonal mean plots with box plot RMSE
fig = plt.figure(figsize =(10, 7))
ax1 = fig.add_subplot(121)
ax1.set_xlabel('\u0394RMSE = $RMSE_{mapped}$ - $RMSE_{unmapped}$')
ax1.set_ylabel('lattitude [\u00b0]')
ax1.plot(rmse_lon_dif.values,rmse_lon_dif.lat.values)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_title("RMSE zonal mean")
ax1.set_title('a.',fontweight = 'bold', loc='left')

ax2 = fig.add_subplot(122)
plt.boxplot(data_rmse,  showmeans=True, showfliers = False)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')
ax2.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax2.set_ylabel('RMSE')
ax2.set_title("RMSE box plot")
ax2.set_title('b.',fontweight = 'bold', loc='left')
plt.savefig('../plot1/global/rmse_lon_bp_'+model+'.png')

#zonal mean plots with box plot KGE
fig= plt.figure(figsize =(10, 7))
ax1 = fig.add_subplot(121)
ax1.set_xlabel('\u0394KGE = $KGE_{mapped}$ - $KGE_{unmapped}$')
ax1.set_ylabel('lattitude [\u00b0]')
ax1.plot(kge_lon_dif.values,kge_lon_dif.lat.values)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_title("KGE zonal mean")
ax1.set_title('a.',fontweight = 'bold', loc='left')

ax2 = fig.add_subplot(122)
plt.boxplot(data_kge, showmeans=True, showfliers = False)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')
ax2.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax2.set_ylabel('KGE')
ax2.set_title('KGE box plot')
ax2.set_title('b.',fontweight = 'bold', loc='left')
plt.savefig('../plot1/global/kge_lon_bp_'+model+'.png')


#zonal mean plots with box plot RMSE seasonal
fig = plt.figure(figsize =(10,7))
ax = fig.add_subplot(111)    # The big subplot
ax.set_xticks([])
ax.set_yticks([])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax.set_xlabel('\u0394RMSE = $RMSE_{mapped}$ - $RMSE_{unmapped}$', labelpad = 35)
ax.set_ylabel('lattitude [\u00b0]', labelpad = 35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax1.set_title( str(obsseas.season.isel(season = 0).values))
ax1.set_title('a.',fontweight = 'bold', loc='left')
ax1.plot(lon_seas[0].values,lon_seas[0].lat.values)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')


ax2.set_title( str(obsseas.season.isel(season = 1).values))
ax2.set_title('b.',fontweight = 'bold', loc='left')
ax2.plot(lon_seas[1].values,lon_seas[1].lat.values)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')

ax3.set_title(str(obsseas.season.isel(season = 2).values))
ax3.set_title('c.',fontweight = 'bold', loc='left')
ax3.plot(lon_seas[2].values,lon_seas[2].lat.values)
ax3.spines['bottom'].set_color('0.75')
ax3.spines['top'].set_color('0.75') 
ax3.spines['right'].set_color('0.75')
ax3.spines['left'].set_color('0.75')

ax4.set_title( str(obsseas.season.isel(season = 3).values))
ax4.set_title('d.',fontweight = 'bold', loc='left')
ax4.plot(lon_seas[3].values,lon_seas[3].lat.values)
ax4.spines['bottom'].set_color('0.75')
ax4.spines['top'].set_color('0.75') 
ax4.spines['right'].set_color('0.75')
ax4.spines['left'].set_color('0.75')

plt.tight_layout()

plt.savefig('../plot1/global/rmse_lon_season_'+model+'.png')

#zonal mean plots with box plot KGE seasonal
fig = plt.figure(figsize =(10,7))
ax = fig.add_subplot(111)    # The big subplot
ax.set_xticks([])
ax.set_yticks([])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax.set_xlabel('\u0394KGE = $KGE_{mapped}$ - $KGE_{unmapped}$',labelpad = 35)
ax.set_ylabel('lattitude [\u00b0]', labelpad = 35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax1.set_title( str(obsseas.season.isel(season = 0).values))
ax1.set_title('a.',fontweight = 'bold', loc='left')
ax1.plot(lon_kg_seas[0].values,lon_kg_seas[0].lat.values)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')


ax2.set_title( str(obsseas.season.isel(season = 1).values))
ax2.set_title('b.',fontweight = 'bold', loc='left')
ax2.plot(lon_kg_seas[1].values,lon_kg_seas[1].lat.values)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')

ax3.set_title(str(obsseas.season.isel(season = 2).values))
ax3.set_title('c.',fontweight = 'bold', loc='left')
ax3.plot(lon_kg_seas[2].values,lon_kg_seas[2].lat.values)
ax3.spines['bottom'].set_color('0.75')
ax3.spines['top'].set_color('0.75') 
ax3.spines['right'].set_color('0.75')
ax3.spines['left'].set_color('0.75')

ax4.set_title( str(obsseas.season.isel(season = 3).values))
ax4.set_title('d.',fontweight = 'bold', loc='left')
ax4.plot(lon_kg_seas[3].values,lon_kg_seas[3].lat.values)
ax4.spines['bottom'].set_color('0.75')
ax4.spines['top'].set_color('0.75') 
ax4.spines['right'].set_color('0.75')
ax4.spines['left'].set_color('0.75')

plt.tight_layout()

plt.savefig('../plot1/global/kge_lon_season_'+model+'.png')


rmse_umpd_seas= xr.open_dataset('../output/rmse_umpd_global_seas_'+model+'.nc')
rmse_mpd_seas= xr.open_dataset('../output/rmse_mpd_global_seas_'+model+'.nc')
rmse_umpdcor_seas= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model+'.nc')
rmse_mpdcor_seas= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model+'.nc')

kge_umpd_seas= xr.open_dataset('../output/kge_umpd_global_seas_'+model + '.nc')
kge_mpd_seas= xr.open_dataset('../output/kge_mpd_global_seas_'+model + '.nc')
kge_umpdcor_seas= xr.open_dataset('../output/kge_umpdcor_global_seas_'+model + '.nc')
kge_mpdcor_seas= xr.open_dataset('../output/kge_mpdcor_global_seas_'+model + '.nc')



fig = plt.figure(figsize =(10, 7))



for i in range(4):
    alp = ['a.', 'b.', 'c.', 'd.']
    pos = [[0.06,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.06,0.06,0.4,0.4],[0.55,0.06,0.4,0.4]]
    
    kge1 = kge_umpd_seas.kge.isel(season =i).stack(z = ("lat","lon"))
    kge2 = kge_mpd_seas.kge.isel(season =i).stack(z = ("lat","lon"))
    kge3 = kge_umpdcor_seas.kge.isel(season =i).stack(z = ("lat","lon"))
    kge4 = kge_mpdcor_seas.kge.isel(season =i).stack(z = ("lat","lon"))

    kge1 = kge1.dropna(dim = 'z')
    kge2 = kge2.dropna(dim = 'z')
    kge3 = kge3.dropna(dim = 'z')
    kge4 = kge4.dropna(dim = 'z')
    
    kge_seas1 = [kge1,kge2, kge3,kge4]
    
    
    
    ax1 = fig.add_axes(pos[i])
    plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
    ax1.spines['bottom'].set_color('0.75')
    ax1.spines['top'].set_color('0.75') 
    ax1.spines['right'].set_color('0.75')
    ax1.spines['left'].set_color('0.75')
    ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
    ax1.set_ylabel('KGE')
    ax1.set_title( str(kge_umpd_seas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
plt.savefig('../plot1/global/kge_bp_season_'+model+'.png') 

fig = plt.figure(figsize =(10, 7))

lon_kg_seas = []

for i in range(4):
    alp = ['a.', 'b.', 'c.', 'd.']
    pos = [[0.06,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.06,0.06,0.4,0.4],[0.55,0.06,0.4,0.4]]
    
    rmse1 = rmse_umpd_seas.rmse.isel(season =i).stack(z = ("lat","lon"))
    rmse2 = rmse_mpd_seas.rmse.isel(season =i).stack(z = ("lat","lon"))
    rmse3 = rmse_umpdcor_seas.rmse.isel(season =i).stack(z = ("lat","lon"))
    rmse4 = rmse_mpdcor_seas.rmse.isel(season =i).stack(z = ("lat","lon"))
    
    rmse1 = rmse1.dropna(dim = 'z')
    rmse2 = rmse2.dropna(dim = 'z')
    rmse3 = rmse3.dropna(dim = 'z')
    rmse4 = rmse4.dropna(dim = 'z')
    
    rmse_seas1 = [rmse1,rmse2, rmse3, rmse4]
    
    ax1 = fig.add_axes(pos[i])
    plt.boxplot(rmse_seas1,  showmeans=True,showfliers = False)
    ax1.spines['bottom'].set_color('0.75')
    ax1.spines['top'].set_color('0.75') 
    ax1.spines['right'].set_color('0.75')
    ax1.spines['left'].set_color('0.75')
    ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
    ax1.set_ylabel('RMSE')
    ax1.set_title( str(rmse_umpd_seas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
plt.savefig('../plot1/global/rmse_bp_season_'+model+'.png') 

    




