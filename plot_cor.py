#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl

model = 'cesm2'
realisation = 'r1i1p1f1'

rmse_umpd = xr.open_dataset('../output/rmse_umpd_global_'+model+'.nc')
rmse_umpdcor = xr.open_dataset('../output/rmse_umpdcor_global_'+model+'.nc')
rmse_mpd = xr.open_dataset('../output/rmse_mpd_global_'+model+'.nc')
rmse_mpdcor = xr.open_dataset('../output/rmse_mpdcor_global_'+model+'.nc')

rmse_mpd = rmse_mpd.sel(lat = slice(84,-56))
rmse_umpd = rmse_umpd.sel(lat = slice(84,-56))

rmse_difumpd  = rmse_umpdcor.rmse - rmse_umpd.rmse
rmse_difmpd  = rmse_mpdcor.rmse - rmse_mpd.rmse

kge_umpd = xr.open_dataset('../output/kge_umpd_global_'+model+'.nc')
kge_umpdcor = xr.open_dataset('../output/kge_umpdcor_global_'+model+'.nc')
kge_mpd = xr.open_dataset('../output/kge_mpd_global_'+model+'.nc')
kge_mpdcor = xr.open_dataset('../output/kge_mpdcor_global_'+model+'.nc')

kge_delta = kge_umpdcor.kge - kge_umpd.kge
kge_deltampd = kge_mpdcor.kge - kge_mpd.kge

kge_delta = kge_delta.sel(lat = slice(84,-56))
kge_deltampd = kge_deltampd.sel(lat = slice(84,-56))

kge_delta = kge_delta.reindex(lat=list(reversed(kge_delta.lat)))
kge_deltampd = kge_deltampd.reindex(lat=list(reversed(kge_deltampd.lat)))

kge_delta= kge_delta.where(kge_delta != 0)
kge_deltampd= kge_deltampd.where(kge_deltampd != 0)

rmse_umpd_seas = xr.open_dataset('../output/rmse_umpd_global_seas_'+model + '.nc')
rmse_mpd_seas = xr.open_dataset('../output/rmse_mpd_global_seas_'+model + '.nc')
rmse_umpdcor_seas = xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model + '.nc')
rmse_mpdcor_seas = xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model + '.nc')

rmse_dif_seas = rmse_mpdcor_seas.rmse - rmse_mpd_seas.rmse
rmse_dif_seas = rmse_dif_seas.sel(lat = slice(84,-56))
rmse_dif_seas = rmse_dif_seas.reindex(lat=list(reversed(rmse_dif_seas.lat)))

rmse_dif_seas = rmse_dif_seas.where(rmse_dif_seas != 0)

kge_umpd_seas = xr.open_dataset('../output/kge_umpd_global_seas_'+model + '.nc')
kge_mpd_seas = xr.open_dataset('../output/kge_mpd_global_seas_'+model + '.nc')
kge_umpdcor_seas = xr.open_dataset('../output/kge_umpdcor_global_seas_'+model + '.nc')
kge_mpdcor_seas = xr.open_dataset('../output/kge_mpdcor_global_seas_'+model + '.nc')

kge_dif_seas = kge_mpdcor_seas.kge - kge_mpd_seas.kge
kge_dif_seas = kge_dif_seas.sel(lat = slice(84,-56))
kge_dif_seas = kge_dif_seas.reindex(lat=list(reversed(kge_dif_seas.lat)))

kge_dif_seas = kge_dif_seas.where(kge_dif_seas != 0)

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

rmse_difumpd = rmse_difumpd.reindex(lat=list(reversed(rmse_difumpd.lat)))
rmse_difmpd = rmse_difmpd.reindex(lat=list(reversed(rmse_difmpd.lat)))

rmse_difumpd = rmse_difumpd.where(rmse_difumpd !=0)
rmse_difmpd = rmse_difmpd.where(rmse_difmpd !=0)
    
lats = rmse_difmpd.lat.values
lons = rmse_difmpd.lon.values
lons, lats = np.meshgrid(lons,lats)

fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_axes([0.05,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,rmse_difumpd.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_title('Root Mean Squared Error (not mapped)')
ax1.set_title('a.',fontweight = 'bold', loc='left')


cbax = fig.add_axes([0.05, 0.2, 0.45, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{corrected}$ - $RMSE_{not-corrected}$', fontsize=9)
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
im1 = m.pcolormesh(lons,lats,kge_delta.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax2.set_title('Kling-Gupta Efficiency (not mapped)')
ax2.set_title('b.',fontweight = 'bold', loc='left')


cbax2 = fig.add_axes([0.53, 0.2, 0.45, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{corrected}$ - $KGE_{not-corrected}$', fontsize=9)
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

plt.savefig('../plot1/global/metrics_difcor_nomap_'+model+'.png')



fig = plt.figure(figsize =(10, 7))

ax1 = fig.add_axes([0.05,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,rmse_difmpd.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_title('Root Mean Squared Error')
ax1.set_title('a.',fontweight = 'bold', loc='left')


cbax = fig.add_axes([0.05, 0.2, 0.45, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{corrected}$ - $RMSE_{not-corrected}$', fontsize=9)
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

ax2 = fig.add_axes([0.53,0.2,0.45,0.5])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax2
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,kge_deltampd.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax2.set_title('Kling-Gupta Efficiency')
ax2.set_title('b.',fontweight = 'bold', loc='left')

cbax2 = fig.add_axes([0.53, 0.2, 0.45, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{corrected}$ - $KGE_{not-corrected}$', fontsize=9)
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)


cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax2.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax2.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax2.transAxes)
plt.savefig('../plot1/global/metrics_difcor_map_'+model+'.png')


lon_seas = []

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
    
    rmse_dif_si = rmse_dif_seas.isel(season = i)
    
    a = rmse_dif_si.mean('lon')
    lon_seas.append(a)
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution=None)
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons,lats,rmse_dif_si.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
    ax1.set_title( str(obsseas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax = fig.add_axes([0.25, 0.1, 0.5, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{corrected}$ - $RMSE_{not-corrected}$', fontsize=9)
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

plt.savefig('../plot1/global/rmse_difcor_map_season_'+model+'.png') 


    


fig = plt.figure(figsize =(10, 7))
ax = fig.add_axes([0.05,0.05,0.95,0.95])
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
    
    kge_dif_si = kge_dif_seas.isel(season = i)
    
    a = kge_dif_si.mean('lon')
    lon_kg_seas.append(a)
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution=None)
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons,lats,kge_dif_si.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
    ax1.set_title( str(kge_dif_seas.season.isel(season = i).values))
    ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax2 = fig.add_axes([0.25, 0.1, 0.5, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{corrected}$ - $KGE_{not-corrected}$', fontsize=9)
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
plt.savefig('../plot1/global/kge_difcor_map_season_'+model+'.png') 







