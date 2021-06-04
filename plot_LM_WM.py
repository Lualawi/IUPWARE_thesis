#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl

model2 = 'cesm2'
model1 = 'ukesm'
realisation = 'r1i1p1f1'



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
cmap.set_bad(color='0.8')

# colorbar args
values1 = [-0.15,-0.125,-0.10,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.10,0.125,0.15]
values2 = [-0.015,-0.0125,-0.010,-0.0075,-0.005,-0.0025,0,0.0025,0.005,0.0075,0.010,0.0125,0.015]
tick_locs1 = [-1.5,-1.0,-0.5,0,0.5,1.0,1.5]
tick_labels1 = ['-0.15','-0.125','-0.10','-0.075','-0.05','-0.025','0','0.025','0.05','0.075','0.10','0.125','0.15']
tick_labels2 = ['-0.015','-0.0125','-0.010','-0.0075','-0.005','-0.0025','0','0.0025','0.005','0.0075','0.010','0.0125','0.015']
norm1 = mpl.colors.BoundaryNorm(values1,cmap.N)
norm2 = mpl.colors.BoundaryNorm(values2,cmap2.N)

rmse_dif1 = xr.open_dataset('../output/rmse_dif_laf_'+model1+'.nc')
rmse_difcor1 = xr.open_dataset('../output/rmse_difcor_laf_'+model1+'.nc')
kge_delta1 = xr.open_dataset('../output/kge_dif_laf_'+model1 + '.nc')
kge_deltacor1 = xr.open_dataset('../output/kge_difcor_laf_'+model1 + '.nc')
rmse_dif1 = rmse_dif1.reindex(lat=list(reversed(rmse_dif1.lat)))
rmse_difcor1 = rmse_difcor1.reindex(lat=list(reversed(rmse_difcor1.lat)))
kge_delta1 = kge_delta1.reindex(lat=list(reversed(kge_delta1.lat)))
kge_deltacor1 = kge_deltacor1.reindex(lat=list(reversed(kge_deltacor1.lat)))
rmse_dif2 = xr.open_dataset('../output/rmse_dif_laf_'+model2+'.nc')
rmse_difcor2 = xr.open_dataset('../output/rmse_difcor_laf_'+model2+'.nc')
kge_delta2 = xr.open_dataset('../output/kge_dif_laf_'+model2 + '.nc')
kge_deltacor2 = xr.open_dataset('../output/kge_difcor_laf_'+model2 + '.nc')
rmse_dif2 = rmse_dif2.reindex(lat=list(reversed(rmse_dif2.lat)))
rmse_difcor2 = rmse_difcor2.reindex(lat=list(reversed(rmse_difcor2.lat)))
kge_delta2 = kge_delta2.reindex(lat=list(reversed(kge_delta2.lat)))
kge_deltacor2 = kge_deltacor2.reindex(lat=list(reversed(kge_deltacor2.lat)))

lats = rmse_dif1.lat.values
lons = rmse_dif1.lon.values
lons, lats = np.meshgrid(lons,lats)

lats1 = rmse_difcor1.lat.values
lons1 = rmse_difcor1.lon.values
lons1, lats1 = np.meshgrid(lons1,lats1)
pos1 = [[0.1,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.1,0.08,0.4,0.4],[0.55,0.08,0.4,0.4]]
pos2 = [[0.25,0.55,0.5,0.015],[0.25,0.065,0.5,0.015]]

fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_axes(pos1[0])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,rmse_dif1.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_ylabel('RMSE (not corrected)')
ax1.set_title('UKESM1',fontweight = 'bold')
ax1.set_title('a.',fontweight = 'bold', loc='left')

cbax = fig.add_axes(pos2[0])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{WM}$ - $RMSE_{LM}$', fontsize=9)
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
ax2 = fig.add_axes(pos1[1])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax2
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,rmse_dif2.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax2.set_title('CESM2',fontweight = 'bold')
ax2.set_title('b.',fontweight = 'bold', loc='left')




ax3 = fig.add_axes(pos1[2])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax3
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,kge_delta1.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax3.set_ylabel('KGE (not corrected)')
ax3.set_title('c.',fontweight = 'bold', loc='left')

cbax3 = fig.add_axes(pos2[1])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax3,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{WM}$ - $KGE_{LM}$', fontsize=9)
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)
#plot arrows
bluelabel1 = 'Higher efficiency/Higher Skill'
redlabel1 = 'Lower efficiency/Lower Skill'

cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax3.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax3.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax3.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax3.transAxes)



ax4 = fig.add_axes(pos1[3])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax4
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons,lats,kge_delta2.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax4.set_title('d.',fontweight = 'bold', loc='left')


plt.savefig('../plot_jun/metrics_diflaf_nocor_both.png')

fig = plt.figure(figsize =(10, 7))


ax1 = fig.add_axes(pos1[0])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax1
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,rmse_difcor1.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax1.set_ylabel('RMSE')
ax1.set_title('UKESM1',fontweight = 'bold')
ax1.set_title('a.',fontweight = 'bold', loc='left')


cbax = fig.add_axes(pos2[0])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb.ax.set_title('\u0394RMSE = $RMSE_{WM}$ - $RMSE_{LM}$', fontsize=9)
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
ax2 = fig.add_axes(pos1[1])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax2
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,rmse_difcor2.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
ax2.set_title('CESM2',fontweight = 'bold')
ax2.set_title('b.',fontweight = 'bold', loc='left')


ax3 = fig.add_axes(pos1[2])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax3
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,kge_deltacor1.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax3.set_ylabel('KGE')
ax3.set_title('c.',fontweight = 'bold', loc='left')

cbax3 = fig.add_axes(pos2[1])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax3,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE = $KGE_{WM}$ - $KGE_{LM}$', fontsize=9)
cb1.ax.set_xticklabels(tick_labels2)
cb1.ax.tick_params(labelsize=7)
#plot arrows
bluelabel1 = 'Higher efficiency/Higher Skill'
redlabel1 = 'Lower efficiency/Lower Skill'

cb1.ax.text(0.7, -3, bluelabel1,size=8, ha='center', va='center',transform  = cbax3.transAxes)
cb1.ax.text(0.3, -3, redlabel1, size=8, ha='center', va='center',transform  = cbax3.transAxes)

cb1.ax.arrow(0.5, -3.5, 0.4, 0, width=0.2, linewidth=0.1, label=bluelabel1,\
          shape='left', head_width=0.5, head_length=0.03,\
          facecolor=cmap40, edgecolor='k', clip_on=False, transform  = cbax3.transAxes)
cb1.ax.arrow(0.5, -3.5, -0.4, 0, width=0.2, linewidth=0.1, label=redlabel1,\
          shape='right', head_width=0.5, head_length=0.03,\
          facecolor=cmap_40, edgecolor='k', clip_on=False, transform  = cbax3.transAxes)



ax4 = fig.add_axes(pos1[3])
m = Basemap(projection='kav7',lon_0=0,resolution=None)
m.ax = ax4
m.drawmapboundary(color='0.8', fill_color='0.8')
im1 = m.pcolormesh(lons1,lats1,kge_deltacor2.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
ax4.set_title('d.',fontweight = 'bold', loc='left')

plt.savefig('../plot_jun/metrics_diflaf_cor_both.png')


alp = ['a.', 'b.', 'c.', 'd.','e.', 'f.', 'g.', 'h.']
pos1 = [[0.1,0.76,0.4,0.2],[0.55,0.76,0.4,0.2],[0.1,0.54,0.4,0.2],[0.55,0.54,0.4,0.2], [0.1,0.32,0.4,0.2],[0.55,0.32,0.4,0.2],[0.1,0.1,0.4,0.2],[0.55,0.1,0.4,0.2]]

rmse_dif_seas1 = xr.open_dataset('../output/rmse_difcor_laf_seas_'+model1+'.nc')
rmse_dif_seas1 = rmse_dif_seas1.reindex(lat=list(reversed(rmse_dif_seas1.lat)))

kge_dif_seas1 = xr.open_dataset('../output/kge_difcor_laf_seas_'+model1+'.nc')
kge_dif_seas1 = kge_dif_seas1.reindex(lat=list(reversed(kge_dif_seas1.lat)))

rmse_dif_seas2 = xr.open_dataset('../output/rmse_difcor_laf_seas_'+model2+'.nc')
rmse_dif_seas2 = rmse_dif_seas2.reindex(lat=list(reversed(rmse_dif_seas2.lat)))

kge_dif_seas2 = xr.open_dataset('../output/kge_difcor_laf_seas_'+model2+'.nc')
kge_dif_seas2 = kge_dif_seas2.reindex(lat=list(reversed(kge_dif_seas2.lat)))


fig = plt.figure(figsize =(10, 14))
ax = fig.add_axes([0.05,0.05,0.9,0.9])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

lon_kg_seas = []

for i in range(8):
  if i%2 == 0:
 	  kge_dif_si = kge_dif_seas1.isel(season = int(i/2))
  if i%2 !=0:
	  kge_dif_si = kge_dif_seas2.isel(season = int((i-1)/2))

  ax1 = fig.add_axes(pos1[i])
  m = Basemap(projection='kav7',lon_0=0,resolution=None)
  m.ax = ax1
  m.drawmapboundary(color='0.8', fill_color='0.8')
  im1 = m.pcolormesh(lons1,lats1,kge_dif_si.kge.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
  ylabel = 'KGE - '+ str(kge_dif_si.season.values)
  if i == 0:
	  ax1.set_title('UKESM1',fontweight = 'bold')
	  ax1.set_ylabel(ylabel)
  elif i==1:
	  ax1.set_title('CESM2',fontweight = 'bold')
  elif i%2 == 0:
	  ax1.set_ylabel(ylabel)
	
  ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax2 = fig.add_axes([0.25, 0.08, 0.5, 0.015])
cb1 = mpl.colorbar.ColorbarBase(ax=cbax2,cmap=cmap2,
                           norm=norm2,
                           ticks = values2,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')
cb1.ax.set_title('\u0394KGE  = $KGE_{WM}$ - $KGE_{LM}$')
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

plt.savefig('../plot_jun/kge_diflaf_cor_season_both.png') 

fig = plt.figure(figsize =(10, 14))
ax = fig.add_axes([0.05,0.05,0.9,0.9])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

lon_kg_seas = []

for i in range(8):
  if i%2 == 0:
 	  rmse_dif_si = rmse_dif_seas1.isel(season = int(i/2))
  elif i%2 !=0:
	  rmse_dif_si = rmse_dif_seas2.isel(season = int((i-1)/2))
    
  ax1 = fig.add_axes(pos1[i])
  m = Basemap(projection='kav7',lon_0=0,resolution=None)
  m.ax = ax1
  m.drawmapboundary(color='0.8', fill_color='0.8')
  im1 = m.pcolormesh(lons1,lats1,rmse_dif_si.rmse.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
  ylabel = 'RMSE - '+ str(rmse_dif_si.season.values)
  if i == 0:
	  ax1.set_title('UKESM1',fontweight = 'bold')
	  ax1.set_ylabel(ylabel)
  elif i==1:
	  ax1.set_title('CESM2',fontweight = 'bold')
  elif i%2 == 0:
	  ax1.set_ylabel(ylabel)
	
  ax1.set_title(alp[i],fontweight = 'bold', loc='left')
    
cbax = fig.add_axes([0.25, 0.08, 0.5, 0.015])
cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                           norm=norm1,
                           ticks = values1,
                           spacing='proportional',
                           orientation='horizontal',
                           extend='both')


cb.ax.set_title('\u0394RMSE = $RMSE_{WM}$ - $RMSE_{LM}$', fontsize=9)
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

plt.savefig('../plot_jun/rmse_diflaf_cor_season_both.png') 