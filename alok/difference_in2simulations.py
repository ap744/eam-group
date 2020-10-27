# This code compare GEOS-Chem model and DEFRA sites ammonia 
# Please contact Alok Pandey ap744@leicester.ac.uk for any further clarifications or details

#import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from sklearn.preprocessing import StandardScaler
import datetime
import xarray as xr
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.cm as cm
import glob
from scipy import stats
#from bootstrap import rma
from scipy.stats import gaussian_kde

#Use todays date to plot files - good practice to save files name with date
Today_date=datetime.datetime.now().strftime("%Y%m%d")

###Different cmap options
# cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')
# cmap = cm.jet
cmap = cm.rainbow
#cmap = cm.YlOrRd

path_naei_run = "/scratch/uptrop/ap744/compare_2_GC_run/geosfp_eu_naei_iccw_SpeciesConcJul2016.nc4"
path_scaleNH3_run = "/scratch/uptrop/ap744/compare_2_GC_run/geosfp_eu_scale_nh3_emis_speciesConcJul2016.nc4"
path_scaleNH3_minus_naei = "/scratch/uptrop/ap744/compare_2_GC_run/speciesconc_scaleNH3_minus_naei_iccw.nc4"

os.chdir('/home/a/ap744/scratch_alok/shapefiles/GBP_shapefile')
Europe_shape = r'GBR_adm1.shp'
Europe_map = ShapelyFeature(Reader(Europe_shape).geometries(),
                               ccrs.PlateCarree(), edgecolor='black',facecolor='none')
print ('Shapefile_read')

naei_run = xr.open_dataset(path_naei_run)
#print (naei_run)
GC_ammonia_naei = naei_run.SpeciesConc_NH3
print (GC_ammonia_naei)
GC_surface_ammonia_naei = GC_ammonia_naei.isel(time=0,lev=0)
print (GC_surface_ammonia_naei)

scaleNH3_run = xr.open_dataset(path_scaleNH3_run)
#print (scaleNH3_run)
GC_ammonia_scaleNH3_run = scaleNH3_run.SpeciesConc_NH3
print (GC_ammonia_scaleNH3_run)
GC_surface_ammonia_scaleNH3_run = GC_ammonia_scaleNH3_run.isel(time=0,lev=0)
print (GC_surface_ammonia_scaleNH3_run)

scaleNH3_minus_naei = xr.open_dataset(path_scaleNH3_minus_naei)
#print (scaleNH3_minus_naei)
GC_ammonia_scaleNH3_minus_naei = scaleNH3_minus_naei.SpeciesConc_NH3
print (GC_ammonia_scaleNH3_minus_naei)
GC_surface_ammonia_scaleNH3_minus_naei = GC_ammonia_scaleNH3_minus_naei.isel(time=0,lev=0)
print (GC_surface_ammonia_scaleNH3_minus_naei)


fig1 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
# plt.suptitle(title_list, fontsize = 35, y=0.96)

ax = plt.subplot(231);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-15, 40, 32, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonia_naei.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1e-9,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'naei_iccw',
											'orientation':'horizontal'})

ax.set_title('naei_iccw',fontsize=25)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
ax.tick_params(labelsize='large')
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'difference_ncfiles_naei_iccw.png',bbox_inches='tight')
plt.show()

fig2 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
# plt.suptitle(title_list, fontsize = 35, y=0.96)

ax = plt.subplot(231);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-15, 40, 32, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonia_scaleNH3_run.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1e-9,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'scaleNH3',
											'orientation':'horizontal'})

ax.set_title('scaleNH3',fontsize=25)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
ax.tick_params(labelsize='large')
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'difference_ncfiles_scaleNH3.png',bbox_inches='tight')
plt.show()


fig3 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
# plt.suptitle(title_list, fontsize = 35, y=0.96)

ax = plt.subplot(231);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-15, 40, 32, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonia_scaleNH3_minus_naei.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1e-9,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'scaleNH3_minus_naei',
											'orientation':'horizontal'})

ax.set_title('scaleNH3_minus_naei',fontsize=25)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
ax.tick_params(labelsize='large')
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'difference_ncfiles_scaleNH3_minus_naei.png',bbox_inches='tight')
plt.show()

