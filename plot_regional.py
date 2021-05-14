#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regionmask
import cartopy.crs as ccrs
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap

model = 'cesm2'
realisation = 'r1i1p1f1'

era = xr.open_dataset('../CMIP/era_mask1.nc')
era_mask = era.t2m.isel(time =474)

obs = xr.open_dataset('../CMIP/Tair_merge.nc')
obseas = obs.groupby('time.season').mean('time')

rmses_dif_clim = pd.read_csv('../output/rmse_dif_regional_clim_'+model+'.csv')
rmses_dif_jja= pd.read_csv('../output/rmse_dif_regional_jja_'+model+'.csv')
rmses_dif_son= pd.read_csv('../output/rmse_dif_regional_son_'+model+'.csv')
rmses_dif_djf= pd.read_csv('../output/rmse_dif_regional_djf_'+model+'.csv')
rmses_dif_mam= pd.read_csv('../output/mse_dif_regional_mam_'+model+'.csv')
rmses_difcor_clim= pd.read_csv('../output/rmse_difcor_regional_clim_'+model+'.csv')

kges_dif_jja= pd.read_csv('../output/kge_dif_regional_jja_'+model+'.csv')
kges_dif_son= pd.read_csv('../output/kge_dif_regional_son_'+model+'.csv')
kges_dif_djf= pd.read_csv('../output/kge_dif_regional_djf_'+model+'.csv')
kges_dif_mam= pd.read_csv('../output/kge_dif_regional_mam_'+model+'.csv')

rmses_difcor_jja= pd.read_csv('../output/rmse_difcor_regional_jja_'+model+'.csv')
rmses_difcor_son= pd.read_csv('../output/rmse_difcor_regional_son_'+model+'.csv')
rmses_difcor_djf= pd.read_csv('../output/rmse_difcor_regional_djf_'+model+'.csv')
rmses_difcor_mam= pd.read_csv('../output/mse_difcor_regional_mam_'+model+'.csv')

kges_difcor_jja= pd.read_csv('../output/kge_difcor_regional_jja_'+model+'.csv')
kges_difcor_son= pd.read_csv('../output/kge_difcor_regional_son_'+model+'.csv')
kges_difcor_djf= pd.read_csv('../output/kge_difcor_regional_djf_'+model+'.csv')
kges_difcor_mam= pd.read_csv('../output/kge_difcor_regional_mam_'+model+'.csv')

kge_dif= pd.read_csv('../output/kge_dif_regional_'+model+'.csv')
kge_difcor= pd.read_csv('../output/kge_difcor_regional_'+model+'.csv')

data = [rmses_dif_clim, rmses_difcor_clim]
data1 = [ kge_dif, kge_difcor]
data2 = [rmses_difcor_jja, rmses_difcor_son, rmses_difcor_mam, rmses_difcor_djf]
data3 = [kges_difcor_jja, kges_difcor_son, kges_difcor_mam, kges_difcor_djf]

lon = np.arange(-179.5, 180)
lat = np.arange(-89.5, 90)
mask = regionmask.defined_regions.ar6.land.mask(lon,lat)
mask =  mask.stack(z = ("lat","lon"))

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
values1 = [-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5]
values2 = [-0.15,-0.125,-0.10,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.10,0.125,0.15]
tick_locs1 = [-1.5,-1.0,-0.5,0,0.5,1.0,1.5]
tick_labels1 = ['-1.5','-1.25','-1.0','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1.0','1.25','1.5']
tick_labels2 = ['-0.15','-0.125','-0.10','-0.075','-0.05','-0.025','0','0.025','0.05','0.075','0.10','0.125','0.15']
norm1 = mpl.colors.BoundaryNorm(values1,cmap.N)
norm2 = mpl.colors.BoundaryNorm(values2,cmap2.N)


lats = mask.lat.values
lons = mask.lon.values
lons, lats = np.meshgrid(lons,lats)

for i in range(len(data)):
    index = np.arange(1,44)
    rmse1 = pd.DataFrame(data = data[i], index = index)
    rmse2 = pd.DataFrame(data = data1[i], index = index)

    print (rmse1.values[0])

    mask1 = mask.copy(deep = True)
    mask1 = mask1.where(mask1 <43)
    
    mask2 = mask.copy(deep = True)
    mask2 = mask1.where(mask1 <43)

    for j in range(len(mask1)):
        if (~np.isnan(mask1[j])):
            if (int(mask1[j].values) == 0):
                mask1[j] = np.nan
            else:
                mask1[j] = float(rmse1.values[int(mask1[j].values)])
        if (~np.isnan(mask2[j])):
            if (int(mask2[j].values) == 0):
                mask2[j] = np.nan
            else:
                mask2[j] = float(rmse2.values[int(mask2[j].values)])
       
    mask1_unst = mask1.unstack("z")
    mask2_unst = mask1.unstack("z")

    
    #using random file to shape coastlines
    mask1_unst = mask1_unst.where(era_mask > 0)
    mask2_unst = mask2_unst.where(era_mask > 0)
    
    lats = mask1_unst.lat.values
    lons = mask1_unst.lon.values
    lons, lats = np.meshgrid(lons,lats)
    
    lats1 = mask2_unst.lat.values
    lons1 = mask2_unst.lon.values
    lons1, lats1 = np.meshgrid(lons1,lats1)
    
    if (i==0):
#       ax.set_title('unmapped-mapped RMSE clim')
        fig= plt.figure(figsize = (10,7))
        
        ax1 = fig.add_axes([0.05,0.05,0.45,0.95])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax1
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons1,lats1,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax1.set_title('Root Mean Squared Error ')
        ax1.set_title('a.',fontweight = 'bold', loc='left')
        
        cbax = fig.add_axes([0.05, 0.27, 0.45, 0.015])
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

        
        ax2 = fig.add_axes([0.53,0.05,0.45,0.95])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax2
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons1,lats1,mask2_unst.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
        ax2.set_title('Kling-Gupta Efficiency')
        ax2.set_title('b.',fontweight = 'bold', loc='left')
        
        cbax2 = fig.add_axes([0.53, 0.27, 0.45, 0.015])
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
        

        plt.savefig('../plot1/regional/metrics_difmap_cor_'+model+'.png')

    else:

        fig= plt.figure(figsize = (10,7))
        
        ax1 = fig.add_axes([0.05,0.05,0.45,0.95])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax1
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax1.set_title('Root Mean Squared Error (not corrected)')
        ax1.set_title('a.',fontweight = 'bold', loc='left')
        cbax = fig.add_axes([0.05, 0.27, 0.45, 0.015])
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
        
        ax2 = fig.add_axes([0.53,0.05,0.45,0.95])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax2
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask2_unst.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
        ax2.set_title('Kling-Gupta Efficiency (not corrected)')
        ax2.set_title('b.',fontweight = 'bold', loc='left')
        
        cbax2 = fig.add_axes([0.53, 0.27, 0.45, 0.015])
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

        
        
        plt.savefig('../plot1/regional/metrics_difmap_nocor_'+model+'.png')


    
fig = plt.figure(figsize =(10,7))
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
    
    index = np.arange(1,44)
    rmse1 = pd.DataFrame(data = data2[i], index = index)
    
    mask1 = mask.copy(deep = True)
    mask1 = mask1.where(mask1 <43)

    for j in range(len(mask1)-1):
        if (~np.isnan(mask1[j])):
            if (int(mask1[j].values) == 0):
                mask1[j] = np.nan
            else:
                mask1[j] = float(rmse1.values[int(mask1[j].values)])
       
    mask1_unst = mask1.unstack("z")
    mask1_unst = mask1_unst.where(era_mask > 0)  
    
    lats = mask1_unst.lat.values
    lons = mask1_unst.lon.values
    lons, lats = np.meshgrid(lons,lats)
    
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution='c')
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons1,lats1,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
    ax1.set_title( str(obseas.isel(season=i).season.values))
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
plt.savefig('../plot1/regional/rmse_difmap_cor_season_'+model+'.png')

fig = plt.figure(figsize =(10,7))
ax = fig.add_axes([0.05,0.05,0.95,0.95])
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

for i in range(4):  
    alp = ['a.', 'b.', 'c.', 'd.']
    pos = [[0.1,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.1,0.13,0.4,0.4],[0.55,0.13,0.4,0.4]]
    
    index = np.arange(1,44)
    kge1 = pd.DataFrame(data = data3[i], index = index)
    
    mask1 = mask.copy(deep = True)
    mask1 = mask1.where(mask1 <43)

    for j in range(len(mask1)-1):
        if (~np.isnan(mask1[j])):
            if (int(mask1[j].values) == 0):
                mask1[j] = np.nan
            else:
                mask1[j] = float(kge1.values[int(mask1[j].values)])
       
    mask1_unst = mask1.unstack("z")
    mask1_unst = mask1_unst.where(era_mask > 0)  
    
    lats = mask1_unst.lat.values
    lons = mask1_unst.lon.values
    lons, lats = np.meshgrid(lons,lats)
    
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution='c')
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons1,lats1,mask1_unst.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
    ax1.set_title( str(obseas.isel(season=i).season.values))
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
plt.savefig('../plot1/regional/kge_difmap_cor_season_'+model+'.png')
