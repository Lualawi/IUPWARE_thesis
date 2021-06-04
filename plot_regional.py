# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 14:08:12 2021

@author: Luwi
"""

#! usr/bin/env python 

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regionmask
import cartopy.crs as ccrs
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap

model1 = 'ukesm'
model2 = 'cesm2'
realisation = 'r1i1p1f2'

era = xr.open_dataset('../CMIP/era_mask1.nc')
era_mask = era.t2m.isel(time =474)

obs = xr.open_dataset('../CMIP/Tair_merge.nc')
obseas = obs.groupby('time.season').mean('time')

rmses_dif_clim1 = pd.read_csv('../output/rmse_dif_regional_clim_'+model1+'.csv')
rmses_dif_jja1= pd.read_csv('../output/rmse_dif_regional_jja_'+model1+'.csv')
rmses_dif_son1= pd.read_csv('../output/rmse_dif_regional_son_'+model1+'.csv')
rmses_dif_djf1= pd.read_csv('../output/rmse_dif_regional_djf_'+model1+'.csv')
rmses_dif_mam1= pd.read_csv('../output/mse_dif_regional_mam_'+model1+'.csv')
rmses_difcor_clim1= pd.read_csv('../output/rmse_difcor_regional_clim_'+model1+'.csv')

kges_dif_jja1= pd.read_csv('../output/kge_dif_regional_jja_'+model1+'.csv')
kges_dif_son1= pd.read_csv('../output/kge_dif_regional_son_'+model1+'.csv')
kges_dif_djf1= pd.read_csv('../output/kge_dif_regional_djf_'+model1+'.csv')
kges_dif_mam1= pd.read_csv('../output/kge_dif_regional_mam_'+model1+'.csv')

rmses_difcor_jja1= pd.read_csv('../output/rmse_difcor_regional_jja_'+model1+'.csv')
rmses_difcor_son1= pd.read_csv('../output/rmse_difcor_regional_son_'+model1+'.csv')
rmses_difcor_djf1= pd.read_csv('../output/rmse_difcor_regional_djf_'+model1+'.csv')
rmses_difcor_mam1= pd.read_csv('../output/mse_difcor_regional_mam_'+model1+'.csv')

kges_difcor_jja1= pd.read_csv('../output/kge_difcor_regional_jja_'+model1+'.csv')
kges_difcor_son1= pd.read_csv('../output/kge_difcor_regional_son_'+model1+'.csv')
kges_difcor_djf1= pd.read_csv('../output/kge_difcor_regional_djf_'+model1+'.csv')
kges_difcor_mam1= pd.read_csv('../output/kge_difcor_regional_mam_'+model1+'.csv')

kge_dif1= pd.read_csv('../output/kge_dif_regional_'+model1+'.csv')
kge_difcor1= pd.read_csv('../output/kge_difcor_regional_'+model1+'.csv')

rmses_dif_clim2 = pd.read_csv('../output/rmse_dif_regional_clim_'+model2+'.csv')
rmses_dif_jja2= pd.read_csv('../output/rmse_dif_regional_jja_'+model2+'.csv')
rmses_dif_son2= pd.read_csv('../output/rmse_dif_regional_son_'+model2+'.csv')
rmses_dif_djf2= pd.read_csv('../output/rmse_dif_regional_djf_'+model2+'.csv')
rmses_dif_mam2= pd.read_csv('../output/mse_dif_regional_mam_'+model2+'.csv')
rmses_difcor_clim2= pd.read_csv('../output/rmse_difcor_regional_clim_'+model2+'.csv')

kges_dif_jja2= pd.read_csv('../output/kge_dif_regional_jja_'+model2+'.csv')
kges_dif_son2= pd.read_csv('../output/kge_dif_regional_son_'+model2+'.csv')
kges_dif_djf2= pd.read_csv('../output/kge_dif_regional_djf_'+model2+'.csv')
kges_dif_mam2= pd.read_csv('../output/kge_dif_regional_mam_'+model2+'.csv')

rmses_difcor_jja2= pd.read_csv('../output/rmse_difcor_regional_jja_'+model2+'.csv')
rmses_difcor_son2= pd.read_csv('../output/rmse_difcor_regional_son_'+model2+'.csv')
rmses_difcor_djf2= pd.read_csv('../output/rmse_difcor_regional_djf_'+model2+'.csv')
rmses_difcor_mam2= pd.read_csv('../output/mse_difcor_regional_mam_'+model2+'.csv')

kges_difcor_jja2= pd.read_csv('../output/kge_difcor_regional_jja_'+model2+'.csv')
kges_difcor_son2= pd.read_csv('../output/kge_difcor_regional_son_'+model2+'.csv')
kges_difcor_djf2= pd.read_csv('../output/kge_difcor_regional_djf_'+model2+'.csv')
kges_difcor_mam2= pd.read_csv('../output/kge_difcor_regional_mam_'+model2+'.csv')

kge_dif2= pd.read_csv('../output/kge_dif_regional_'+model2+'.csv')
kge_difcor2= pd.read_csv('../output/kge_difcor_regional_'+model2+'.csv')


rs1 = [rmses_dif_clim1, rmses_difcor_clim1]
rs2 = [rmses_dif_clim2,rmses_difcor_clim1]
kg1 = [ kge_dif1, kge_difcor1]
kg2 = [ kge_dif2, kge_difcor2]
data2 = [rmses_difcor_djf1,rmses_difcor_djf2, rmses_difcor_jja1,rmses_difcor_jja2, rmses_difcor_mam1,rmses_difcor_mam2, rmses_difcor_son1,rmses_difcor_son2]
data3 = [kges_difcor_djf1,kges_difcor_djf2, kges_difcor_jja1,kges_difcor_jja2, kges_difcor_mam1,kges_difcor_mam2, kges_difcor_son1,kges_difcor_son2]

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

pos1 = [[0.1,0.57,0.4,0.4],[0.55,0.57,0.4,0.4],[0.1,0.08,0.4,0.4],[0.55,0.08,0.4,0.4]]
pos2 = [[0.25,0.55,0.5,0.015],[0.25,0.065,0.5,0.015]]


for i in range(len(rs1)):
    index = np.arange(1,44)
    rmse1 = pd.DataFrame(data = rs1[i], index = index)
    rmse2 = pd.DataFrame(data = rs2[i], index = index)
    kge1 = pd.DataFrame(data = kg1[i], index = index)
    kge2 = pd.DataFrame(data = kg2[i], index = index)


    mask1 = mask.copy(deep = True)
    mask1 = mask1.where(mask1 <43)
    
    mask2 = mask.copy(deep = True)
    mask2 = mask2.where(mask1 <43)
    
    mask3 = mask.copy(deep = True)
    mask3 = mask3.where(mask1 <43)
    
    mask4 = mask.copy(deep = True)
    mask4 = mask4.where(mask1 <43)

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
        if (~np.isnan(mask3[j])):
            if (int(mask3[j].values) == 0):
                mask3[j] = np.nan
            else:
                mask3[j] = float(kge1.values[int(mask3[j].values)])
        if (~np.isnan(mask4[j])):
            if (int(mask4[j].values) == 0):
                mask4[j] = np.nan
            else:
                mask4[j] = float(kge2.values[int(mask4[j].values)])
       
    mask1_unst = mask1.unstack("z")
    mask2_unst = mask2.unstack("z")
    mask3_unst = mask3.unstack("z")
    mask4_unst = mask4.unstack("z")

    
    #using random file to shape coastlines
    mask1_unst = mask1_unst.where(era_mask > 0)
    mask2_unst = mask2_unst.where(era_mask > 0)
    mask3_unst = mask3_unst.where(era_mask > 0)
    mask4_unst = mask4_unst.where(era_mask > 0)
    
    lats = mask1_unst.lat.values
    lons = mask1_unst.lon.values
    lons, lats = np.meshgrid(lons,lats)
    
    
    if (i==0):
#       ax.set_title('unmapped-mapped RMSE clim')
        fig= plt.figure(figsize = (10,7))
        
        ax1 = fig.add_axes(pos1[0])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax1
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
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

        
        ax2 = fig.add_axes(pos1[1])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax2
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask2_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax2.set_title('CESM2',fontweight = 'bold')
        ax2.set_title('b.',fontweight = 'bold', loc='left')
        
        
        ax3 = fig.add_axes(pos1[2])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax3
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask3_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax3.set_ylabel('KGE (not corrected)')
        ax3.set_title('c.',fontweight = 'bold', loc='left')

        cbax2 = fig.add_axes(pos2[1])
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
        
        ax4 = fig.add_axes(pos1[3])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax4
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask4_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax4.set_title('d.',fontweight = 'bold', loc='left')


        plt.savefig('../plot_jun/metrics_difmap_nocor_reg.png')

    else:

        fig= plt.figure(figsize = (10,7))
        
        ax1 = fig.add_axes(pos1[0])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax1
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
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

        
        ax2 = fig.add_axes(pos1[1])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax2
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask2_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax2.set_title('CESM2',fontweight = 'bold')
        ax2.set_title('b.',fontweight = 'bold', loc='left')
        
        
        ax3 = fig.add_axes(pos1[2])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax3
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask3_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax3.set_ylabel('KGE')
        ax3.set_title('c.',fontweight = 'bold', loc='left')

        cbax2 = fig.add_axes(pos2[1])
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
        
        ax4 = fig.add_axes(pos1[3])
        m = Basemap(projection='kav7',lon_0=0,resolution='c')
        m.ax = ax4
        m.drawmapboundary(color='0.8', fill_color='0.8')
        im1 = m.pcolormesh(lons,lats,mask4_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
        ax4.set_title('d.',fontweight = 'bold', loc='left')
        
        
        plt.savefig('../plot_jun/metrics_difmap_cor_reg.png')


    
fig = plt.figure(figsize =(10,14))


alp = ['a.', 'b.', 'c.', 'd.','e.', 'f.', 'g.', 'h.']
pos = [[0.1,0.76,0.4,0.2],[0.55,0.76,0.4,0.2],[0.1,0.54,0.4,0.2],[0.55,0.54,0.4,0.2], [0.1,0.32,0.4,0.2],[0.55,0.32,0.4,0.2],[0.1,0.1,0.4,0.2],[0.55,0.1,0.4,0.2]]

for i in range(8):  
    
    if i%2 == 0:
 	    seas = obseas.isel(season = int(i/2)).season.values
    elif i%2 !=0:
	    seas = obseas.isel(season = int((i-1)/2)).season.values
    
    print(str(seas))
    index = np.arange(1,44)
    rmse = pd.DataFrame(data = data2[i], index = index)
    
    mask1 = mask.copy(deep = True)
    mask1 = mask1.where(mask1 <43)

    for j in range(len(mask1)-1):
        if (~np.isnan(mask1[j])):
            if (int(mask1[j].values) == 0):
                mask1[j] = np.nan
            else:
                mask1[j] = float(rmse.values[int(mask1[j].values)])
       
    mask1_unst = mask1.unstack("z")
    mask1_unst = mask1_unst.where(era_mask > 0)  
    
    lats = mask1_unst.lat.values
    lons = mask1_unst.lon.values
    lons, lats = np.meshgrid(lons,lats)
    
    
    ax1 = fig.add_axes(pos[i])
    m = Basemap(projection='kav7',lon_0=0,resolution='c')
    m.ax = ax1
    m.drawmapboundary(color='0.8', fill_color='0.8')
    im1 = m.pcolormesh(lons,lats,mask1_unst.values,shading='flat',cmap=cmap,norm = norm1,latlon=True)
    ylabel = 'RMSE - '+ str(seas)
    
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
plt.savefig('../plot_jun/rmse_difmap_cor_season_reg.png')

fig = plt.figure(figsize =(10,14))


for i in range(8):  
    if i%2 == 0:
 	    seas = obseas.isel(season = int(i/2)).season.values
    elif i%2 !=0:
	    seas = obseas.isel(season = int((i-1)/2)).season.values

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
    im1 = m.pcolormesh(lons,lats,mask1_unst.values,shading='flat',cmap=cmap2,norm = norm2,latlon=True)
    ylabel = 'KGE - '+ str(seas)
    
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
plt.savefig('../plot_jun/kge_difmap_cor_season_reg.png')
