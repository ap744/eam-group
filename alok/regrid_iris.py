import os
import sys
import iris
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import pandas as pd
import math
import calendar
import datetime
#import seaborn as sns # changes some plot settings (lineweights, fonts)
import warnings
import traceback
import os
import netCDF4 as nc4
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
import cartopy.feature as cfeature
from mpl_toolkits.basemap import Basemap, shiftgrid
import cartopy.crs as ccrs
import matplotlib.cm as cm
import cartopy
#cmap = cm.jet
cmap = cm.rainbow
cmap1 = cm.coolwarm
import datetime
from scipy import stats
from scipy.stats import gaussian_kde
Today_date=datetime.datetime.now().strftime("%Y%m%d")

def discrete_cmap(N, base_cmap=None):
	"""Create an N-bin discrete colormap from the specified input map"""
	# Note that if base_cmap is a string or None, you can simply do
	#    return plt.cm.get_cmap(base_cmap, N)
	# The following works for string, None, or a colormap instance:
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0, 1, N))
	cmap_name = base.name + str(N)
	return base.from_list(cmap_name, color_list, N)

def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
		lons and lats are 1-d array while data is 2-D array
		colorbar_min,colorbar_max specifies the minimum and maximum value you want to show with the hottest and coldest color respectively
		tb_lef and tb_bot specifies if you want to have axis labels, True fro yes
	output : a spatial map of the data
	"""
	lons[lons>180]-=360; 
	# lon_b = np.min(lons); lon_e = np.max(lons)
	# lon_b = -180; lon_e = 180
	# lon_b = 65; lon_e = 100
	lon_b = -9; lon_e = 3 #lonW,lonE,
	
	# lat_b = np.min(lats); lat_e = np.max(lats)	
	# lat_b = -90; lat_e = 90
	# lat_b = 5; lat_e = 38
	lat_b = 49; lat_e = 61 #latS,latN
	
	lon_bin = 4; lat_bin = 4
	
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)

	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=24)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=24)
	# Add Coastlines, States, and Country Boundaries
	# map.drawcoastlines(); map.drawcountries() #map.drawstates(); # draw border lines
	map.readshapefile('/scratch/uptrop/ap744/shapefiles/Shapfiles_india/World_shp/World','World',linewidth=2) # to add a shapefile on a map
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	#masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)   #  use 20 color bins, this can be changed
	cmap.set_bad([1,1,1],alpha = 1.0);
	if bad_data:
		#cmap.set_under('w');cmap.set_over('k')
		cmap.set_under('k');cmap.set_over('w')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	return colormesh	


cubes = iris.load('/scratch/uptrop/ap744/python_work/20200724emission_monthly_GC.nc')
#print(cubes)

Jan = cubes[0]
cube1 = Jan
print (Jan.shape)
lat_o = Jan.coord('latitude')
print(lat_o.shape)
lat_o =lat_o.points[:]
print (lat_o,'points')

lon_o = Jan.coord('longitude')
print(lon_o)

lon_o =lon_o.points[:]

nh3_o_mar = Jan.data[4,:,:]
print (nh3_o_mar.shape)


lat,lon = cube1.coord('latitude'),cube1.coord('longitude')
lat_min,lon_min = lat.points.min(),lon.points.min()
lat_max,lon_max = lat.points.max(),lon.points.max()
lat25 = np.arange(lat_min, lat_max, 0.1)
lon25 = np.arange(lon_min, lon_max, 0.1)


result = cube1.interpolate([('latitude', lat25), ('longitude', lon25)],
                           iris.analysis.Linear())

print(result.shape)

lat_n = result.coord('latitude')
print(lat_n.shape)

lon_n = result.coord('longitude')
print(lon_n)

nh3_n_mar = result.data[4,:,:]
print (nh3_n_mar.shape)

#lat_o=np.arange(32.75,61.50,0.25)
#lon_o= np.arange(-15,40.3125,0.3125)


lat_n=np.arange(32.75,61.25,0.1)
lon_n= np.arange(-15,40,0.1)

title_list = 'NH3 Emission regridding comparison'					###########
title_list1 = 'Original'
title_list2 = 'Regridded'

fig = plt.figure(facecolor='White',figsize=[11,7]);pad= 1.1; 
plt.suptitle(title_list, fontsize = 25, y=1.001)

ax = plt.subplot(1,2,1);
plt.title(title_list1, fontsize = 30, y=1)
colormap=cmap; colorbar_min=1e-07;colorbar_max=1e-02## change this accordingly
colormesh_1 = spatial_figure(ax,nh3_o_mar,lon_o,lat_o,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
cbar_ax = fig.add_axes([0.10, 0.02, 0.35, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'Original',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(1,2,2);
plt.title(title_list2, fontsize = 30, y=1)
colormap=cmap; colorbar_min=1e-07;colorbar_max=1e-02## change this accordingly
colormesh_1 = spatial_figure(ax,nh3_n_mar,lon_n,lat_n,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	


cbar_ax = fig.add_axes([0.60, 0.02, 0.35, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'Regrid',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	

plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'regrid_iris_may.png',bbox_inches='tight')   ##########	
plt.show()


