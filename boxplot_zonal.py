#! usr/bin/env python 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

cmap_whole = plt.cm.get_cmap('RdBu_r')
cmap40 = cmap_whole(0.15)

model1 = 'ukesm'
model2 = 'cesm2'

rmse_umpd_seas1= xr.open_dataset('../output/rmse_umpd_global_seas_'+model1+'.nc')
rmse_mpd_seas1= xr.open_dataset('../output/rmse_mpd_global_seas_'+model1+'.nc')
rmse_umpdcor_seas1= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model1+'.nc')
rmse_mpdcor_seas1= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model1+'.nc')

kge_umpd_seas1= xr.open_dataset('../output/kge_umpd_global_seas_'+model1 + '.nc')
kge_mpd_seas1= xr.open_dataset('../output/kge_mpd_global_seas_'+model1 + '.nc')
kge_umpdcor_seas1= xr.open_dataset('../output/kge_umpdcor_global_seas_'+model1 + '.nc')
kge_mpdcor_seas1= xr.open_dataset('../output/kge_mpdcor_global_seas_'+model1 + '.nc')

rmse_umpd_seas2= xr.open_dataset('../output/rmse_umpd_global_seas_'+model2+'.nc')
rmse_mpd_seas2= xr.open_dataset('../output/rmse_mpd_global_seas_'+model2+'.nc')
rmse_umpdcor_seas2= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model2+'.nc')
rmse_mpdcor_seas2= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model2+'.nc')

kge_umpd_seas2= xr.open_dataset('../output/kge_umpd_global_seas_'+model2 + '.nc')
kge_mpd_seas2= xr.open_dataset('../output/kge_mpd_global_seas_'+model2 + '.nc')
kge_umpdcor_seas2= xr.open_dataset('../output/kge_umpdcor_global_seas_'+model2 + '.nc')
kge_mpdcor_seas2= xr.open_dataset('../output/kge_mpdcor_global_seas_'+model2 + '.nc')


alp = ['a.', 'b.', 'c.', 'd.','e.', 'f.', 'g.', 'h.']
pos = [[0.1,0.76,0.4,0.2],[0.55,0.76,0.4,0.2],[0.1,0.52,0.4,0.2],[0.55,0.52,0.4,0.2], [0.1,0.28,0.4,0.2],[0.55,0.28,0.4,0.2],[0.1,0.04,0.4,0.2],[0.55,0.04,0.4,0.2]]

fig = plt.figure(figsize =(10, 14))



kge1 = kge_umpd_seas1.kge.isel(season =0).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas1.kge.isel(season =0).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas1.kge.isel(season =0).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas1.kge.isel(season =0).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[0])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('KGE - '+str(kge_umpd_seas1.season.isel(season = 0).values))
ax1.set_title('UKESM1',fontweight = 'bold')
ax1.set_title(alp[0],fontweight = 'bold', loc='left')
ax1.arrow(1.2, 0.6, 0.6, 0.01, width=0.04, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.05, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 0.6, 0.6, 0.01, width=0.04, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.05, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

kge1 = kge_umpd_seas2.kge.isel(season =0).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas2.kge.isel(season =0).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas2.kge.isel(season =0).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas2.kge.isel(season =0).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[1])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title('CESM2',fontweight = 'bold')
ax1.set_title(alp[1],fontweight = 'bold', loc='left')


kge1 = kge_umpd_seas1.kge.isel(season =1).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas1.kge.isel(season =1).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas1.kge.isel(season =1).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas1.kge.isel(season =1).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[2])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('KGE - '+ str(kge_umpd_seas1.season.isel(season = 1).values))
ax1.set_title(alp[2],fontweight = 'bold', loc='left')
#ax1.arrow(1.2, 0.4, 0.6, 0.01, width=0.01, linewidth=0.1, capstyle = 'butt',\
#          shape='right', head_width=0.02, head_length=0.1,\
#          facecolor=cmap40, edgecolor='k', clip_on=False)
          
#ax1.arrow(3.2, 0.4, 0.6, 0.01, width=0.01, linewidth=0.1,capstyle = 'butt',\
#          shape='right', head_width=0.02, head_length=0.1,\
#          facecolor=cmap40, edgecolor='k', clip_on=False)

kge1 = kge_umpd_seas2.kge.isel(season =1).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas2.kge.isel(season =1).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas2.kge.isel(season =1).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas2.kge.isel(season =1).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[3])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[3],fontweight = 'bold', loc='left')



kge1 = kge_umpd_seas1.kge.isel(season =2).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas1.kge.isel(season =2).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas1.kge.isel(season =2).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas1.kge.isel(season =2).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[4])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('KGE - '+ str(kge_umpd_seas1.season.isel(season = 2).values))
ax1.set_title(alp[4],fontweight = 'bold', loc='left')

ax1.arrow(1.2, 0.9, 0.6, 0.01, width=0.016, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.026, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 0.9, 0.6, 0.01, width=0.016, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.026, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

kge1 = kge_umpd_seas2.kge.isel(season =2).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas2.kge.isel(season =2).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas2.kge.isel(season =2).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas2.kge.isel(season =2).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[5])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[5],fontweight = 'bold', loc='left')

ax1.arrow(1.2, 0.9, 0.6, 0.01, width=0.016, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.026, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 0.9, 0.6, 0.01, width=0.016, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.026, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

kge1 = kge_umpd_seas1.kge.isel(season =3).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas1.kge.isel(season =3).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas1.kge.isel(season =3).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas1.kge.isel(season =3).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[6])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('KGE - '+ str(kge_umpd_seas1.season.isel(season = 3).values))
ax1.set_title(alp[6],fontweight = 'bold', loc='left')

kge1 = kge_umpd_seas2.kge.isel(season =3).stack(z = ("lat","lon"))
kge2 = kge_mpd_seas2.kge.isel(season =3).stack(z = ("lat","lon"))
kge3 = kge_umpdcor_seas2.kge.isel(season =3).stack(z = ("lat","lon"))
kge4 = kge_mpdcor_seas2.kge.isel(season =3).stack(z = ("lat","lon"))

kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

kge_seas1 = [kge1,kge2, kge3,kge4]

ax1 = fig.add_axes(pos[7])
plt.boxplot(kge_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[7],fontweight = 'bold', loc='left')

ax1.arrow(1.2, 0.4, 0.6, 0.01, width=0.02, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.03, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 0.4, 0.6, 0.01, width=0.02, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.03, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
plt.savefig('../plot_jun/kge_bp_seasonal_both.png') 

rmse_umpd_seas1= xr.open_dataset('../output/rmse_umpd_global_seas_'+model1+'.nc')
rmse_mpd_seas1= xr.open_dataset('../output/rmse_mpd_global_seas_'+model1+'.nc')
rmse_umpdcor_seas1= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model1+'.nc')
rmse_mpdcor_seas1= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model1+'.nc')

rmse_umpd_seas1= xr.open_dataset('../output/rmse_umpd_global_seas_'+model1 + '.nc')
rmse_mpd_seas1= xr.open_dataset('../output/rmse_mpd_global_seas_'+model1 + '.nc')
rmse_umpdcor_seas1= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model1 + '.nc')
rmse_mpdcor_seas1= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model1 + '.nc')

rmse_umpd_seas2= xr.open_dataset('../output/rmse_umpd_global_seas_'+model2+'.nc')
rmse_mpd_seas2= xr.open_dataset('../output/rmse_mpd_global_seas_'+model2+'.nc')
rmse_umpdcor_seas2= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model2+'.nc')
rmse_mpdcor_seas2= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model2+'.nc')

rmse_umpd_seas2= xr.open_dataset('../output/rmse_umpd_global_seas_'+model2 + '.nc')
rmse_mpd_seas2= xr.open_dataset('../output/rmse_mpd_global_seas_'+model2 + '.nc')
rmse_umpdcor_seas2= xr.open_dataset('../output/rmse_umpdcor_global_seas_'+model2 + '.nc')
rmse_mpdcor_seas2= xr.open_dataset('../output/rmse_mpdcor_global_seas_'+model2 + '.nc')


fig = plt.figure(figsize =(10, 14))

rmse1 = rmse_umpd_seas1.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas1.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas1.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas1.rmse.isel(season =0).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[0])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('RMSE - '+str(rmse_umpd_seas1.season.isel(season = 0).values))
ax1.set_title('UKESM1',fontweight = 'bold')
ax1.set_title(alp[0],fontweight = 'bold', loc='left')


rmse1 = rmse_umpd_seas2.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas2.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas2.rmse.isel(season =0).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas2.rmse.isel(season =0).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[1])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title('CESM2',fontweight = 'bold')
ax1.set_title(alp[1],fontweight = 'bold', loc='left')



rmse1 = rmse_umpd_seas1.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas1.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas1.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas1.rmse.isel(season =1).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[2])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('RMSE - '+ str(rmse_umpd_seas1.season.isel(season = 1).values))
ax1.set_title(alp[2],fontweight = 'bold', loc='left')
ax1.arrow(1.2, 5, 0.6, -0.1, width=0.1, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 5, 0.6, -0.1, width=0.1, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

rmse1 = rmse_umpd_seas2.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas2.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas2.rmse.isel(season =1).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas2.rmse.isel(season =1).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[3])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[3],fontweight = 'bold', loc='left')
ax1.arrow(1.2, 5, 0.6, -0.1, width=0.1, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 5, 0.6, -0.1, width=0.1, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

rmse1 = rmse_umpd_seas1.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas1.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas1.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas1.rmse.isel(season =2).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[4])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('RMSE - '+ str(rmse_umpd_seas1.season.isel(season = 2).values))
ax1.set_title(alp[4],fontweight = 'bold', loc='left')

ax1.arrow(1.2, 6, 0.6, -0.1, width=0.1, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 6, 0.6, -0.1, width=0.1, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

rmse1 = rmse_umpd_seas2.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas2.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas2.rmse.isel(season =2).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas2.rmse.isel(season =2).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[5])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[5],fontweight = 'bold', loc='left')

ax1.arrow(1.2, 6, 0.6, -0.1, width=0.1, linewidth=0.1, capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)
          
ax1.arrow(3.2, 6, 0.6, -0.1, width=0.1, linewidth=0.1,capstyle = 'butt',\
          shape='right', head_width=0.2, head_length=0.1,\
          facecolor=cmap40, edgecolor='k', clip_on=False)

rmse1 = rmse_umpd_seas1.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas1.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas1.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas1.rmse.isel(season =3).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[6])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_ylabel('RMSE - '+ str(rmse_umpd_seas1.season.isel(season = 3).values))
ax1.set_title(alp[6],fontweight = 'bold', loc='left')

rmse1 = rmse_umpd_seas2.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse2 = rmse_mpd_seas2.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor_seas2.rmse.isel(season =3).stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor_seas2.rmse.isel(season =3).stack(z = ("lat","lon"))

rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

rmse_seas1 = [rmse1,rmse2, rmse3,rmse4]

ax1 = fig.add_axes(pos[7])
plt.boxplot(rmse_seas1,  showmeans=True, showfliers = False)
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax1.set_title(alp[7],fontweight = 'bold', loc='left')
    
plt.savefig('../plot_jun/rmse_bp_seasonal_both.png') 

rmse_dif_cm = xr.open_dataset('../output/rmse_dif_global_'+model1+'.nc')
rmse_difcor_cm = xr.open_dataset('../output/rmse_difcor_global_'+model1+'.nc')

rmse_dif_um = xr.open_dataset('../output/rmse_dif_global_'+model2+'.nc')
rmse_difcor_um = xr.open_dataset('../output/rmse_difcor_global_'+model2+'.nc')

obs = xr.open_dataset('../CMIP/Tair_merge.nc')
obsseas = obs.Tair.groupby('time.season').mean('time')

rmse_umpd1 = xr.open_dataset('../output/rmse_umpd_global_'+model1+'.nc')
rmse_mpd1= xr.open_dataset('../output/rmse_mpd_global_'+model1+'.nc')
rmse_umpdcor1= xr.open_dataset('../output/rmse_umpdcor_global_'+model1+'.nc')
rmse_mpdcor1= xr.open_dataset('../output/rmse_mpdcor_global_'+model1+'.nc')

kge_umpd1= xr.open_dataset('../output/kge_umpd_global_'+model1 + '.nc')
kge_mpd1= xr.open_dataset('../output/kge_mpd_global_'+model1 + '.nc')
kge_umpdcor1= xr.open_dataset('../output/kge_umpdcor_global_'+model1 + '.nc')
kge_mpdcor1= xr.open_dataset('../output/kge_mpdcor_global_'+model1 + '.nc')

rmse1 = rmse_umpd1.rmse.stack(z = ("lat","lon"))
rmse2 = rmse_mpd1.rmse.stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor1.rmse.stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor1.rmse.stack(z = ("lat","lon"))
rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

data_rmse_cm = [rmse1, rmse2, rmse3, rmse4]

#plot box and longitude for KGE
kge1 = kge_umpd1.kge.stack(z = ("lat","lon"))
kge2 = kge_mpd1.kge.stack(z = ("lat","lon"))
kge3 = kge_umpdcor1.kge.stack(z = ("lat","lon"))
kge4 = kge_mpdcor1.kge.stack(z = ("lat","lon"))
kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

data_kge_cm = [kge1, kge2, kge3, kge4]

rmse_lon_umpd1 = rmse_umpd1.rmse.mean(dim = 'lon')
rmse_lon_mpd1 = rmse_mpd1.rmse.mean(dim = 'lon')
rmse_lon_dif1 = (rmse_lon_mpd1 - rmse_lon_umpd1).rename('rmse')

rmse_lon_umpdcor1 = rmse_umpdcor1.rmse.mean(dim = 'lon')
rmse_lon_mpdcor1 = rmse_mpdcor1.rmse.mean(dim = 'lon')
rmse_lon_difcor1 = (rmse_lon_mpdcor1 - rmse_lon_umpdcor1).rename('rmse')

kge_lon_umpd1 = kge_umpd1.kge.mean(dim = 'lon')
kge_lon_mpd1 = kge_mpd1.kge.mean(dim = 'lon')
kge_lon_dif1 = (kge_lon_mpd1 - kge_lon_umpd1).rename('kge')

kge_lon_umpdcor1 = kge_umpdcor1.kge.mean(dim = 'lon')
kge_lon_mpdcor1 = kge_mpdcor1.kge.mean(dim = 'lon')
kge_lon_difcor1 = (kge_lon_mpdcor1 - kge_lon_umpdcor1).rename('kge')


rmse_umpd2 = xr.open_dataset('../output/rmse_umpd_global_'+model2+'.nc')
rmse_mpd2= xr.open_dataset('../output/rmse_mpd_global_'+model2+'.nc')
rmse_umpdcor2= xr.open_dataset('../output/rmse_umpdcor_global_'+model2+'.nc')
rmse_mpdcor2= xr.open_dataset('../output/rmse_mpdcor_global_'+model2+'.nc')

kge_umpd2= xr.open_dataset('../output/kge_umpd_global_'+model2 + '.nc')
kge_mpd2= xr.open_dataset('../output/kge_mpd_global_'+model2 + '.nc')
kge_umpdcor2= xr.open_dataset('../output/kge_umpdcor_global_'+model2 + '.nc')
kge_mpdcor2= xr.open_dataset('../output/kge_mpdcor_global_'+model2 + '.nc')

rmse1 = rmse_umpd2.rmse.stack(z = ("lat","lon"))
rmse2 = rmse_mpd2.rmse.stack(z = ("lat","lon"))
rmse3 = rmse_umpdcor2.rmse.stack(z = ("lat","lon"))
rmse4 = rmse_mpdcor2.rmse.stack(z = ("lat","lon"))
rmse1 = rmse1.dropna(dim = 'z')
rmse2 = rmse2.dropna(dim = 'z')
rmse3 = rmse3.dropna(dim = 'z')
rmse4 = rmse4.dropna(dim = 'z')

data_rmse_um = [rmse1, rmse2, rmse3, rmse4]

#plot box and longitude for KGE
kge1 = kge_umpd2.kge.stack(z = ("lat","lon"))
kge2 = kge_mpd2.kge.stack(z = ("lat","lon"))
kge3 = kge_umpdcor2.kge.stack(z = ("lat","lon"))
kge4 = kge_mpdcor2.kge.stack(z = ("lat","lon"))
kge1 = kge1.dropna(dim = 'z')
kge2 = kge2.dropna(dim = 'z')
kge3 = kge3.dropna(dim = 'z')
kge4 = kge4.dropna(dim = 'z')

data_kge_um = [kge1, kge2, kge3, kge4]

rmse_lon_umpd2 = rmse_umpd2.rmse.mean(dim = 'lon')
rmse_lon_mpd2 = rmse_mpd2.rmse.mean(dim = 'lon')
rmse_lon_dif2 = (rmse_lon_mpd2 - rmse_lon_umpd2).rename('rmse')

rmse_lon_umpdcor2 = rmse_umpdcor2.rmse.mean(dim = 'lon')
rmse_lon_mpdcor2 = rmse_mpdcor2.rmse.mean(dim = 'lon')
rmse_lon_difcor2 = (rmse_lon_mpdcor2 - rmse_lon_umpdcor2).rename('rmse')

kge_lon_umpd2 = kge_umpd2.kge.mean(dim = 'lon')
kge_lon_mpd2 = kge_mpd2.kge.mean(dim = 'lon')
kge_lon_dif2 = (kge_lon_mpd2 - kge_lon_umpd2).rename('kge')

kge_lon_umpdcor2 = kge_umpdcor2.kge.mean(dim = 'lon')
kge_lon_mpdcor2 = kge_mpdcor2.kge.mean(dim = 'lon')
kge_lon_difcor2 = (kge_lon_mpdcor2 - kge_lon_umpdcor2).rename('kge')

alon = [rmse_lon_difcor1,rmse_lon_difcor2, kge_lon_difcor1,kge_lon_difcor2]
alon2 = []
for dat in alon:
  datval = dat.values
  #print('orig:',' '.join(str(x) for x in datval))
  dat2 = []
  for i in range(len(datval)):
    if i <2 or i > len(datval)-3:
      dat2.append(0)
    else:
      avg = np.mean([datval[i-2],datval[i-1],datval[i],datval[i+1],datval[i+2]])
      dat2.append(avg)
  dat3 = 0*dat.copy()
  dat3.values = dat2
  alon2.append(dat3)

rmse_lon_difcor1.to_pandas().to_csv('orig.csv')
rmse_lon_difcor1 = alon2[0]
rmse_lon_difcor1.to_pandas().to_csv('smooth.csv')
rmse_lon_difcor2 = alon2[1]
kge_lon_difcor1 = alon2[2]
kge_lon_difcor2 = alon2[3]
  

#zonal mean plots with box plot RMSE
fig = plt.figure(figsize =(10, 7))
ax1 = fig.add_subplot(121)
ax1.set_xlabel('\u0394RMSE = $RMSE_{mapped}$ - $RMSE_{unmapped}$')
ax1.set_ylabel('lattitude [\u00b0]')
ax1.plot(rmse_lon_difcor1.values,rmse_lon_difcor1.lat.values,color="blue",label = "CESM2")
ax1.plot(rmse_lon_difcor2.values,rmse_lon_difcor2.lat.values,color="red",label = "UKESM1")
ax1.legend()
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_title("RMSE zonal mean")
ax1.set_title('a.',fontweight = 'bold', loc='left')

ax2 = fig.add_subplot(222)
plt.boxplot(data_rmse_cm,  showmeans=True, showfliers = False)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')
ax2.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax2.set_ylabel('RMSE')
ax2.set_title("RMSE box plot (CESM2)")
ax2.set_title('b.',fontweight = 'bold', loc='left')

ax3 = fig.add_subplot(224)
plt.boxplot(data_rmse_um,  showmeans=True, showfliers = False)
ax3.spines['bottom'].set_color('0.75')
ax3.spines['top'].set_color('0.75') 
ax3.spines['right'].set_color('0.75')
ax3.spines['left'].set_color('0.75')
ax3.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax3.set_ylabel('RMSE')
ax3.set_title("RMSE box plot (UKESM1)")
ax3.set_title('c.',fontweight = 'bold', loc='left')

plt.tight_layout()
plt.savefig('../plot_jun/rmse_lon_bp_both.png')

#zonal mean plots with box plot KGE
fig= plt.figure(figsize =(10, 7))
ax1 = fig.add_subplot(121)
ax1.set_xlabel('\u0394KGE = $KGE_{mapped}$ - $KGE_{unmapped}$')
ax1.set_ylabel('lattitude [\u00b0]')
ax1.plot(kge_lon_difcor1.values,kge_lon_difcor1.lat.values,color="blue",label = "CESM2")
ax1.plot(kge_lon_difcor2.values,kge_lon_difcor2.lat.values,color="red",label = "UKESM1")
ax1.legend()
ax1.spines['bottom'].set_color('0.75')
ax1.spines['top'].set_color('0.75') 
ax1.spines['right'].set_color('0.75')
ax1.spines['left'].set_color('0.75')
ax1.set_title("KGE zonal mean")
ax1.set_title('a.',fontweight = 'bold', loc='left')

ax2 = fig.add_subplot(222)
plt.boxplot(data_kge_cm, showmeans=True, showfliers = False)
ax2.spines['bottom'].set_color('0.75')
ax2.spines['top'].set_color('0.75') 
ax2.spines['right'].set_color('0.75')
ax2.spines['left'].set_color('0.75')
ax2.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax2.set_ylabel('KGE')
ax2.set_title('KGE box plot (CESM2)')
ax2.set_title('b.',fontweight = 'bold', loc='left')

ax3 = fig.add_subplot(224)
plt.boxplot(data_kge_um, showmeans=True, showfliers = False)
ax3.spines['bottom'].set_color('0.75')
ax3.spines['top'].set_color('0.75') 
ax3.spines['right'].set_color('0.75')
ax3.spines['left'].set_color('0.75')
ax3.set_xticklabels(['umpd', 'mpd', 'umpdcor', 'mpdcor'])
ax3.set_ylabel('KGE')
ax3.set_title('KGE box plot (UKESM1)')
ax3.set_title('c.',fontweight = 'bold', loc='left')

plt.tight_layout()

plt.savefig('../plot_jun/kge_lon_bp_both.png')



