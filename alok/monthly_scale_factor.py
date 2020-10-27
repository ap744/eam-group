import netCDF4 as nc4
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
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
	# lon_bin = 60; lat_bin = 30
	lon_bin = 4; lat_bin = 4
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)

	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=24)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=24)
	# Add Coastlines, States, and Country Boundaries
	# map.drawcoastlines(); map.drawcountries() #map.drawstates(); # draw border lines, here only coast lines
	map.readshapefile('/scratch/uptrop/ap744/shapefiles/Shapfiles_india/World_shp/World','World', linewidth=2) # to add a shapefile on a map
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	#masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)   #  use 20 color bins, this can be changed
	cmap.set_bad([1,1,1],alpha = 1.0);
	if bad_data:
		cmap.set_under('w');cmap.set_over('k')
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

	#satellite data files
	def create_filename_sat(month,day):
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
		sat_data_files = '/scratch/uptrop/em440/for_Alok/gc_ncdf/satellite_files/ts_08_11.EU.2016'+str(month)+str(day)+'.nc'
		return sat_data_files
	
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
	
	#Satellites daily files 
	sat_data_files = create_filename_sat(month,day)
	
	ncf_sat = nc4.Dataset(sat_data_files,mode='r')
	lat_sat = ncf_sat.variables['LAT'][:]
	lon_sat = ncf_sat.variables['LON'][:]
	nh3_sat = ncf_sat.variables['IJ-AVG-S__NH3'][:]     		#NH3 tracer 'ppbv'
	airdensity_sat = ncf_sat.variables['TIME-SER__AIRDEN'][:]	#Air density 'molecules/cm3'
	bxheight_sat = ncf_sat.variables['BXHGHT-S__BXHEIGHT'][:]	#Grid Box height 'm'
	print (nh3_sat.shape, airdensity_sat.shape, bxheight_sat.shape)
	
	bxheight1_sat = bxheight_sat*100 #Grid Box height 'cm'
	airdensity1_sat = airdensity_sat * bxheight1_sat # Air density 'molecules/cm3' * Grid Box height 'cm' = molecules/cm2
	
	nh3_sat_A = nh3_sat*airdensity1_sat #molucules/cm2
	print (nh3_sat_A.shape)
	
	nh3_sat_B = np.nansum(nh3_sat_A, axis=0) #sum over all model vertical layers
	print (nh3_sat_B.shape)
	
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
	nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)
	print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	
	nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) #kg(NH3)/cm2
	print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	
	nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)           #moles(NO2)/cm2
	print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	
	nh3_emission_total_E = nh3_emission_total_D*NA                      #molecules/cm2
	print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')
	
	return nh3_sat_B,lon_sat,lat_sat, nh3_emission_total_E,lon_emission,lat_emission

def monthly_mean_cal():	
	months=np.arange(1,13)
	time_series = np.empty(len(months));	

	sat_mon_mean = np.empty((len(time_series),115,177))
	sat_mon_mean[:] = np.nan
	
	emission_mon_mean = np.empty((len(time_series),115,177))
	emission_mon_mean[:] = np.nan

	
	for imonth in months:
		sat_nh3_mon_mean = np.empty((115,177))
		sat_nh3_mon_mean[:] = np.nan
		
		emission_nh3_mon_mean = np.empty((115,177))
		emission_nh3_mon_mean[:] = np.nan
				
		days = days_months(imonth)
		for iday in days:	
			nh3_sat_B,lon_sat,lat_sat, nh3_emission_total_E,lon_emission,lat_emission = daily_data(imonth,iday)
			sat_nh3_mon_mean = np.dstack((sat_nh3_mon_mean,nh3_sat_B))
			emission_nh3_mon_mean = np.dstack((emission_nh3_mon_mean,nh3_emission_total_E))
			
			#print (sat_nh3_mon_mean.shape)
		sat_mon_mean[imonth-1,:,:] = np.nanmean(sat_nh3_mon_mean,axis=2)
		emission_mon_mean[imonth-1,:,:] = np.nanmean(emission_nh3_mon_mean,axis=2)
		
	return time_series, lon_sat, lat_sat, sat_mon_mean, emission_mon_mean

time_series, lon_sat, lat_sat, sat_mon_mean, emission_mon_mean = monthly_mean_cal()

print(sat_mon_mean.shape,'!sat_mon_mean')
print(emission_mon_mean.shape,'!emission_mon_mean')


annual_emission = np.nansum(emission_mon_mean, axis=0)
print (annual_emission.shape, 'annual_emission.shape')
print (np.nanmax(annual_emission), 'max',np.nanmin(annual_emission), 'min')


Jan_scale_factor = emission_mon_mean[0,:,:]/annual_emission
Feb_scale_factor = emission_mon_mean[1,:,:]/annual_emission
Mar_scale_factor = emission_mon_mean[2,:,:]/annual_emission
Apr_scale_factor = emission_mon_mean[3,:,:]/annual_emission
May_scale_factor = emission_mon_mean[4,:,:]/annual_emission
Jun_scale_factor = emission_mon_mean[5,:,:]/annual_emission
Jul_scale_factor = emission_mon_mean[6,:,:]/annual_emission
Aug_scale_factor = emission_mon_mean[7,:,:]/annual_emission
Sep_scale_factor = emission_mon_mean[8,:,:]/annual_emission
Oct_scale_factor = emission_mon_mean[9,:,:]/annual_emission
Nov_scale_factor = emission_mon_mean[10,:,:]/annual_emission
Dec_scale_factor = emission_mon_mean[11,:,:]/annual_emission

print (np.nanmax(Jan_scale_factor), 'max',np.nanmin(Jan_scale_factor), 'min')

title_list = 'Monthly scale factor (i_mon_emission/annual_emission)'					###########
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
fig = plt.figure(facecolor='White',figsize=[21,30]);pad= 1.1; 
plt.suptitle(title_list, fontsize = 38, y=1.001)

ax = plt.subplot(4,3,1);
plt.title(title_list1, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2## change this accordingly
colormesh_1 = spatial_figure(ax,Jan_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,2);
plt.title(title_list2, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2## change this accordingly
colormesh_1 = spatial_figure(ax,Feb_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,3);
plt.title(title_list3, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2 ## change this accordingly
colormesh_1 = spatial_figure(ax,Mar_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
ax = plt.subplot(4,3,4);
plt.title(title_list4, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,Apr_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,5);
plt.title(title_list5, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,May_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,6);
plt.title(title_list6, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,Jun_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,7);
plt.title(title_list7, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,Jul_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,8);
plt.title(title_list8, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,Aug_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,9);
plt.title(title_list9, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2  ## change this accordingly
colormesh_1 = spatial_figure(ax,Sep_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,10);
plt.title(title_list10, fontsize = 30, y=1)
colormap=cmap;colorbar_min=0.02;colorbar_max=0.2 ## change this accordingly
colormesh_1 = spatial_figure(ax,Oct_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,11);
plt.title(title_list11, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2 ## change this accordingly
colormesh_1 = spatial_figure(ax,Nov_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(4,3,12);
plt.title(title_list12, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.02;colorbar_max=0.2 ## change this accordingly
colormesh_1 = spatial_figure(ax,Dec_scale_factor,lon_sat,lat_sat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

cbar_ax = fig.add_axes([0.10, 0.02, 0.80, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'scale factor',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'monthly_scale_factor.png',bbox_inches='tight')   ##########	
plt.show()


