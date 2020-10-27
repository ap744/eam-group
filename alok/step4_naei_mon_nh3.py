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
import netCDF4 as nc4
import matplotlib
import cartopy.feature as cfeature
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.cm as cm
import cartopy
#cmap = cm.jet
cmap = cm.rainbow
cmap1 = cm.coolwarm
import datetime
from scipy import stats
from scipy.stats import gaussian_kde
Today_date=datetime.datetime.now().strftime("%Y%m%d")
from iris.coords import DimCoord
from iris.cube import Cube

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
	#lons[lons>180]-=360; 
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

def days_months(months):
	if months in [1,3,5,7,8,10,12]:
		days=np.arange(1,32)
	elif months in [4,6,9,11]:
		days=np.arange(1,31)
	else:
		days=np.arange(1,30) #as feb 2016 has 29 days
	return days

def daily_data(month,day):
	NA=6.022e23   #molecules/mol
	mNH3=17.0      #g(NO2)/mol
	mair=28.97    #g(air)/mol

	#emissions data files
	def create_filename_emission(month,day):
		"""
		# define file names to be read in each loop
		"""
		if month<=9:
			month = '0'+str(month);
		else:
			month =str(month);
		print (month, 'iN @ Month')
		if day<=9:
			day = '0'+str(day);
		else:
			day =str(day);
		print (day, 'iN @ Day')
		emissions_data_files = '/scratch/uptrop/em440/for_Alok/gc_ncdf/emissions/HEMCO_diagnostics.2016'+str(month)+str(day)+'0000.nc'
		return emissions_data_files

	#emissions daily files
	emission_data_files = create_filename_emission(month,day)
	
	ncf_emission = nc4.Dataset(emission_data_files,mode='r')
	lat_emission = ncf_emission.variables['lat'][:]
	lon_emission = ncf_emission.variables['lon'][:]
	area_emission = ncf_emission.variables['AREA'][:]  
	nh3_emission_total = ncf_emission.variables['EmisNH3_Total'][:]  
	nh3_emission_anthro = ncf_emission.variables['EmisNH3_Anthro'][:]
	print (area_emission.shape, 'area_emission.shape')
	print (nh3_emission_total.shape, 'nh3_emission_total.shape')
	
	nh3_emission_total_A = nh3_emission_total[0,:,:,:]
	print (nh3_emission_total_A.shape, 'nh3_emission_total_A.shape')
	#sum over all vertical layers
	nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)			#kg
	print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission) 			#kg(NH3)/m2
	#nh3_emission_total_E = nh3_emission_total_C*1000 						#g/m2/month as calculating fora a month
	
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 	#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)          #moles(NO2)/cm2
	#print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	#nh3_emission_total_E = nh3_emission_total_D*NA                     #molecules/cm2
	#print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')	
	
	return lat_emission, lon_emission, nh3_emission_total_B

def monthly_mean_cal():	
	months=np.arange(1,13)
	time_series = np.empty(len(months));	

	emission_mon_mean = np.empty((len(time_series),115,177))
	emission_mon_mean[:] = np.nan
	
	for imonth in months:

		emission_nh3_mon_mean = np.empty((115,177))
		emission_nh3_mon_mean[:] = np.nan
				
		days = days_months(imonth)
		for iday in days:	
			lat_emission, lon_emission, nh3_emission_total_E = daily_data(imonth,iday)
			emission_nh3_mon_mean = np.dstack((emission_nh3_mon_mean,nh3_emission_total_E))
			#print (emission_nh3_mon_mean.shape)
		emission_mon_mean[imonth-1,:,:] = np.nanmean(emission_nh3_mon_mean,axis=2)
		
	return time_series, lat_emission, lon_emission, emission_mon_mean

time_series, lat, lon, emission_mon_mean = monthly_mean_cal()
print(emission_mon_mean.shape,'!emission_mon_mean')


#regridding using iris
lat_min,lon_min = np.nanmin(lat),np.nanmin(lon)
lat_max,lon_max = np.nanmax(lat),np.nanmax(lon)
lat01 = np.arange(lat_min, lat_max, 0.1)
lon01 = np.arange(lon_min, lon_max, 0.1)

latitude = DimCoord(lat,
					standard_name='latitude',
					units='degrees')
longitude = DimCoord(lon,
					standard_name='longitude',
					units='degrees')
time = DimCoord(np.linspace(1, 12, 12),
					standard_name='time',
					units='month')

cube2 = Cube(emission_mon_mean,
					dim_coords_and_dims=[(latitude, 1),
										(longitude, 2),
										(time, 0)])

regridded_data_emission = cube2.interpolate([('latitude', lat01), ('longitude', lon01)],
                           iris.analysis.Linear())
print(regridded_data_emission.shape, 'regrid_emission shape')

lat_n = regridded_data_emission.coord('latitude')
#print(lat_n.shape)
lat_n =lat_n.points[:]
lon_n = regridded_data_emission.coord('longitude')
#print(lon_n)
lon_n =lon_n.points[:]
regridded_data_emission = regridded_data_emission.data[:]
annual_emission = np.nansum(regridded_data_emission.data[:], axis=0)

lat_n_min,lon_n_min = np.nanmin(lat_n),np.nanmin(lon_n)
lat_n_max,lon_n_max = np.nanmax(lat_n),np.nanmax(lon_n)
print (lat_n_min, 'lat_min_scale_factor')
print (lon_n_min, 'lon_min_scale_factor')
print (lat_n_max, 'lat_max_scale_factor')
print (lon_n_max, 'lon_max_scale_factor')

#monthly scale factor
Jan_scale_factor = regridded_data_emission[0,:,:]/annual_emission #unit - unitless(kg/kg from GC model)
Feb_scale_factor = regridded_data_emission[1,:,:]/annual_emission
Mar_scale_factor = regridded_data_emission[2,:,:]/annual_emission
Apr_scale_factor = regridded_data_emission[3,:,:]/annual_emission
May_scale_factor = regridded_data_emission[4,:,:]/annual_emission
Jun_scale_factor = regridded_data_emission[5,:,:]/annual_emission
Jul_scale_factor = regridded_data_emission[6,:,:]/annual_emission
Aug_scale_factor = regridded_data_emission[7,:,:]/annual_emission
Sep_scale_factor = regridded_data_emission[8,:,:]/annual_emission
Oct_scale_factor = regridded_data_emission[9,:,:]/annual_emission
Nov_scale_factor = regridded_data_emission[10,:,:]/annual_emission
Dec_scale_factor = regridded_data_emission[11,:,:]/annual_emission

print (Jan_scale_factor.shape, 'Jan_scale_factor.shape')

naei_nh3_file = nc4.Dataset('/scratch/uptrop/em440/for_Alok/naei_nh3/NAEI_total_NH3_0.1x0.1_2016.nc',mode='r')
lat_naei = naei_nh3_file.variables['lat'][:]
lon_naei = naei_nh3_file.variables['lon'][:]
naei_nh3 = naei_nh3_file.variables['NH3'][:] 	#unit g/m2/yr
naei_area = naei_nh3_file.variables['area'][:] 	#unit m2

naei_nh3_area = (naei_nh3 * naei_area )/1000 # g/m2/yr * m2 = g/yr --> g/yr/1000 --->kg/yr
naei_nh3_area_mon = naei_nh3_area/12 # kg/month
naei_nh3_area_mon[naei_nh3_area_mon<100] =np.nan
#print (naei_nh3.shape, 'naei_nh3.shape')
lat_naei_min,lon_naei_min = np.nanmin(lat_naei),np.nanmin(lon_naei)
lat_naei_max,lon_naei_max = np.nanmax(lat_naei),np.nanmax(lon_naei)
#print (lat_naei_min, 'lat_min_naei')
#print (lon_naei_min, 'lon_min_naei')
#print (lat_naei_max, 'lat_max_naei')
#print (lon_naei_max, 'lon_max_naei')
#print (naei_nh3.shape, 'naei_nh3.shape')
#print (lat_n, 'lat_scale_factor')
#print (lon_n, 'lon_scale_factor')
#print (lat_naei, 'lat_naei')
#print (lon_naei, 'lon_naei')


print (naei_nh3.shape, 'naei_nh3.shape')

Jan_scale_factor_uk = Jan_scale_factor[165:,46:200]
print (Jan_scale_factor_uk.shape, 'Jan_scale_factor_uk')
naei_nh3_uk = naei_nh3_area_mon[0:120,:]
print (naei_nh3_uk.shape)
lat_uk = lat01[165:]
lon_uk = lon01[46:200]
lat_naei_uk = lat_naei[0:120]
lon_naei_uk = lon_naei[:]
print (lat_uk.shape, lat_uk, 'lat_scale_factor_uk')
print (lat_naei_uk.shape, lat_naei_uk, 'lat_naei_uk')
print (lon_uk.shape, lon_uk, 'lon_scale_factor_uk')
print (lon_naei_uk.shape, lon_naei_uk, 'lon_naei_uk')

Feb_scale_factor_uk = Feb_scale_factor[165:,46:200]
Mar_scale_factor_uk = Mar_scale_factor[165:,46:200]
Apr_scale_factor_uk = Apr_scale_factor[165:,46:200]
May_scale_factor_uk = May_scale_factor[165:,46:200]
Jun_scale_factor_uk = Jun_scale_factor[165:,46:200]
Jul_scale_factor_uk = Jul_scale_factor[165:,46:200]
Aug_scale_factor_uk = Aug_scale_factor[165:,46:200]
Sep_scale_factor_uk = Sep_scale_factor[165:,46:200]
Oct_scale_factor_uk = Oct_scale_factor[165:,46:200]
Nov_scale_factor_uk = Nov_scale_factor[165:,46:200]
Dec_scale_factor_uk = Dec_scale_factor[165:,46:200]

Jan_naei_nh3_emission = Jan_scale_factor_uk *  naei_nh3_uk
Feb_naei_nh3_emission = Feb_scale_factor_uk *  naei_nh3_uk
Mar_naei_nh3_emission = Mar_scale_factor_uk *  naei_nh3_uk
Apr_naei_nh3_emission = Apr_scale_factor_uk *  naei_nh3_uk
May_naei_nh3_emission = May_scale_factor_uk *  naei_nh3_uk
Jun_naei_nh3_emission = Jun_scale_factor_uk *  naei_nh3_uk
Jul_naei_nh3_emission = Jul_scale_factor_uk *  naei_nh3_uk
Aug_naei_nh3_emission = Aug_scale_factor_uk *  naei_nh3_uk
Sep_naei_nh3_emission = Sep_scale_factor_uk *  naei_nh3_uk
Oct_naei_nh3_emission = Oct_scale_factor_uk *  naei_nh3_uk
Nov_naei_nh3_emission = Nov_scale_factor_uk *  naei_nh3_uk
Dec_naei_nh3_emission = Dec_scale_factor_uk *  naei_nh3_uk



print (np.nanmax(Apr_naei_nh3_emission), 'Apr_Max', np.nanmin(Apr_naei_nh3_emission), 'Apr_Min')

title_list = 'NAEI monthly NH$_3$ Emission 0.1 x 0.1'					###########
title_list1 = 'JAN'
title_list2 = 'FEB'
title_list3 = 'MAR'
title_list4 = 'APR'
title_list5 = 'MAY'
title_list6 = 'JUN'
title_list7 = 'JUL'
title_list8 = 'AUG'
title_list9 = 'SEP'
title_list10 = 'OCT'
title_list11 = 'NOV'
title_list12 = 'DEC'

# Plot Figure 
fig = plt.figure(facecolor='White',figsize=[33,25]);pad= 1.1; 
plt.suptitle(title_list, fontsize = 38, y=1.001)

ax = plt.subplot(3,4,1);
plt.title(title_list1, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Jan_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,2);
plt.title(title_list2, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Feb_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,3);
plt.title(title_list3, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Mar_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
ax = plt.subplot(3,4,4);
plt.title(title_list4, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,Apr_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,5);
plt.title(title_list5, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,May_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,6);
plt.title(title_list6, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,Jun_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,7);
plt.title(title_list7, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,Jul_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,8);
plt.title(title_list8, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,Aug_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,9);
plt.title(title_list9, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3 ## change this accordingly
colormesh_1 = spatial_figure(ax,Sep_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,10);
plt.title(title_list10, fontsize = 30, y=1)
colormap=cmap;colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Oct_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,11);
plt.title(title_list11, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Nov_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,12);
plt.title(title_list12, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=1e3## change this accordingly
colormesh_1 = spatial_figure(ax,Dec_naei_nh3_emission,lon_uk,lat_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

cbar_ax = fig.add_axes([0.10, 0.02, 0.80, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'NAEI Emission (kg/mon)',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'step4_naei_monthly_emission11.png',bbox_inches='tight')   ##########	
plt.show()
