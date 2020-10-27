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
	""" Calculate number of days in a month """
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
	area_emission = ncf_emission.variables['AREA'][:]  					#unit m2
	nh3_emission_total = ncf_emission.variables['EmisNH3_Total'][:]  	#Unit kg
	nh3_emission_anthro = ncf_emission.variables['EmisNH3_Anthro'][:] 	#Unit kg
	print (area_emission.shape, 'area_emission.shape')
	print (nh3_emission_total.shape, 'nh3_emission_total.shape')
	
	nh3_emission_total_A = nh3_emission_total[0,:,:,:]
	print (nh3_emission_total_A.shape, 'nh3_emission_total_A.shape')
	#sum over all vertical layers
	nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)			#kg  --> sum over all vertical layers
	print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 		#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)           	#moles(NO2)/cm2
	#print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	#nh3_emission_total_E = nh3_emission_total_D*NA                      	#molecules/cm2
	#print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')
	print (lat_sat.shape, lon_sat.shape, 'lat lon shape of sat file')
	print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')
	return lat_sat,lon_sat,nh3_sat_B, lat_emission,lon_emission,nh3_emission_total_B

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
			lat_sat,lon_sat,nh3_sat_B, lat_emission,lon_emission,nh3_emission_total_B = daily_data(imonth,iday)
			sat_nh3_mon_mean = np.dstack((sat_nh3_mon_mean,nh3_sat_B))
			emission_nh3_mon_mean = np.dstack((emission_nh3_mon_mean,nh3_emission_total_B))
			
			#print (sat_nh3_mon_mean.shape)
		sat_mon_mean[imonth-1,:,:] = np.nanmean(sat_nh3_mon_mean,axis=2)
		emission_mon_mean[imonth-1,:,:] = np.nanmean(emission_nh3_mon_mean,axis=2)
	print (lat_sat.shape, lon_sat.shape, 'lat lon shape of sat file')
	print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')	
	
	return time_series, lat_sat, lon_sat, emission_mon_mean, sat_mon_mean 

time_series, lat, lon, emission_mon_mean, sat_mon_mean = monthly_mean_cal()
print(sat_mon_mean.shape,'!sat_mon_mean.shape')
print(emission_mon_mean.shape,'!emission_mon_mean.shape')

def netcdf4_write(file_name,time_series,lat,lon,emission_nh3_GC,nh3_column_conc_GC):
	"""
	Here is how you would normally create and store data in a netCDF file:
    1) Open/create a netCDF dataset.
    2) Define the dimensions of the data.
    3) Construct netCDF variables using the defined dimensions.
    4) Pass data into the netCDF variables.
    5) Add attributes to the variables and dataset (optional but recommended).
    6) Close the netCDF dataset.
	"""

	# 1) Open/create a netCDF dataset.
	f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	
	"""
	The above line creates a netCDF file called "file_name" in the filepath folder.
	f is a netCDF Dataset object that provides methods for storing data to the file. 
	f also doubles as the root group. A netCDF group is basically a directory or 
	folder within the netCDF dataset. This allows you to organize data as you would 
	in a unix file system. Let's create a group for the heck of it
	"""

	# 2) Define the dimensions of the data.
	"""
	netCDF defines the sizes of all variables in terms of dimensions, so before any 
	variables can be created the dimen sions they use must be created first. A special 
	case, not often used in practice, is that of a scalar variable, which has no dimensions.
	A dimension is created using the Dataset.createDimension method of a Dataset. A Python 
	string is used to set the name of the dimension, and an integer value is used to set 
	the size. To create an unlimited dimension (a dimension that can be appended to), t
	he size value is set to None or 0. In this example, the time dimension is unlimited. 
	In netCDF4 files you can have more than one unlimited dimension, in netCDF 3 files 
	there may be only one, and it must be the first (leftmost) dimension of the variable
	"""
	f.createDimension('time', len(time_series))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	# f.createDimension('height', len(height))

	#3) Construct netCDF variables using the defined dimensions.
	times = f.createVariable('time',np.float64, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	# heights = f.createVariable('height',np.float32, ('height'))
	# ozones = f.createVariable('ozone',np.float32,('time','height','lat','lon'))
	emission_nh3_GCs = f.createVariable('emission_nh3_GC',np.float32,('time','lat','lon'))
	nh3_column_conc_GCs = f.createVariable('nh3_column_conc_GC',np.float32,('time','lat','lon'))
	'''
	prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
	total_precips = f.createVariable('total_precip',np.float32,('time','lat','lon'))
	mean_preps = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
	std_preps = f.createVariable('std_precip',np.float32,('time','lat','lon'))
	MPIs = f.createVariable('MPI',np.float32,('time','lat','lon'))
	'''
	#4) Passing data into variables
	times[:] = time_series
	latitudes[:] = lat
	longitudes[:] = lon
	emission_nh3_GCs[:] = emission_nh3_GC
	nh3_column_conc_GCs[:] = nh3_column_conc_GC
	# heights[:] = height
	# ozones[:] = ozone
	'''
	std_preps[:] = std_prep
	MPIs[:] = MPI
	'''

	# 5) Add attributes to the variables and dataset (optional but recommended).
	"""
	There are two types of attributes in a netCDF file, global and variable. Global attributes provide information
	about the entire dataset, as a whole. Variable attributes provide information about one of the variables. Global
	attributes are set by assigning values to Dataset instance variables. Variable attributes are set by assigning
	values to Variable instances variables. Attributes can be strings, numbers or sequences. 
	"""
	# Global Attributes
	"""
	title       : Whats in the file
	institution : Where it was produced
	source      : How it was produced e.g. model version, instrument type
	history     : Audit trail of processing operations
	references  : Pointers to publications or web documentation
	comment     : Miscellaneous
	"""
	f.description = 'Geos-Chem emission daily to monthly mean datasets'
	# f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alok Pandey at the University of Edinburgh'

	# Variable Attributes
	"""
	units               : mandatory for all variables containing data other than dimensionless numbers. 
	standard_name       : identifies the quantity. Units must be consistent with standard name and any 
						  statistical processing.Standard name does not include coordinate or processing information.
	long_name           : not standardised.
	ancillary variables : a pointer to variables providing metadata about the individual data values e.g. standard error 
						  or data quality information. Numeric data variables may have 
						  ##############################################################################
						  # FillValue, missing_value, valid_max, valid_min, valid_range. missing_value #
						  ##############################################################################
	"""
	times.long_name = 'Month'
	times.units = '1'
	latitudes.long_name = 'latitude'
	latitudes.units = 'degree_north'
	longitudes.long_name = 'longitude'
	longitudes.units = 'degree_east'
	# heights.units = 'meters'

	emission_nh3_GCs.standard_name = 'Geos-Chem NH3 monthly emission'
	emission_nh3_GCs.long_name = 'Geos-Chem NH3 monthly emission'
	emission_nh3_GCs.units = 'kg'
	nh3_column_conc_GCs.standard_name = 'Geos-Chem monthly NH3 column concentration'
	nh3_column_conc_GCs.long_name = 'Geos-Chem monthly NH3 column concentration'
	nh3_column_conc_GCs.units = 'molucules/cm2'
	# ozones.standard_name = 'Ozone concentration'
	# ozones.long_name = 'Ozone hourly concentration'
	# ozones.units = 'ppb'
	'''
	rx5days.standard_name = 'maximum 5-day precip amount'
	rx5days.long_name = 'maximum 5-day precipitation amount'
	rx5days.units = 'mm' 
	'''
	
	f.close()
	print ('nc file write:Done!!!')
	return 1

#write nc file
file_name= Today_date+'monthly_NH3columnConc_and_emission_monthly_kg.nc'
netcdf4_write('/scratch/uptrop/ap744/python_work/'+file_name,time_series, lat, lon, emission_mon_mean, sat_mon_mean)

print (np.nanmax(sat_mon_mean), 'max_sat_mon_mean', np.nanmin(sat_mon_mean), 'min_sat_mon_mean' )
print (np.nanmax(emission_mon_mean), 'max_emission_mon_mean', np.nanmin(emission_mon_mean), 'min_emission_mon_mean' )
#mask low emission mainly over ocean
emission_mon_mean[emission_mon_mean <= 250] = np.nan

#calculation of GC ratio (NH3 column concentration / NH3 emission)
Jan_mon_mean_model_ratio = sat_mon_mean[0,:,:]/emission_mon_mean[0,:,:]
Feb_mon_mean_model_ratio = sat_mon_mean[1,:,:]/emission_mon_mean[1,:,:]
Mar_mon_mean_model_ratio = sat_mon_mean[2,:,:]/emission_mon_mean[2,:,:]
Apr_mon_mean_model_ratio = sat_mon_mean[3,:,:]/emission_mon_mean[3,:,:]
May_mon_mean_model_ratio = sat_mon_mean[4,:,:]/emission_mon_mean[4,:,:]
Jun_mon_mean_model_ratio = sat_mon_mean[5,:,:]/emission_mon_mean[5,:,:]
Jul_mon_mean_model_ratio = sat_mon_mean[6,:,:]/emission_mon_mean[6,:,:]
Aug_mon_mean_model_ratio = sat_mon_mean[7,:,:]/emission_mon_mean[7,:,:]
Sep_mon_mean_model_ratio = sat_mon_mean[8,:,:]/emission_mon_mean[8,:,:]
Oct_mon_mean_model_ratio = sat_mon_mean[9,:,:]/emission_mon_mean[9,:,:]
Nov_mon_mean_model_ratio = sat_mon_mean[10,:,:]/emission_mon_mean[10,:,:]
Dec_mon_mean_model_ratio = sat_mon_mean[11,:,:]/emission_mon_mean[11,:,:]

print (np.nanmax(Jan_mon_mean_model_ratio),np.nanmin(Jan_mon_mean_model_ratio),'max & Min')

title_list = 'Monthly Mean GC ratio (NH$_3$ column conc/emission)'					###########
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
plt.suptitle(title_list, fontsize = 35, y=1.001)

ax = plt.subplot(3,4,1);
plt.title(title_list1, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21## change this accordingly
colormesh_1 = spatial_figure(ax,Jan_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,2);
plt.title(title_list2, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21## change this accordingly
colormesh_1 = spatial_figure(ax,Feb_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,3);
plt.title(title_list3, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21 ## change this accordingly
colormesh_1 = spatial_figure(ax,Mar_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
ax = plt.subplot(3,4,4);
plt.title(title_list4, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,Apr_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,5);
plt.title(title_list5, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,May_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,6);
plt.title(title_list6, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,Jun_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('JJA',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,7);
plt.title(title_list7, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,Jul_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,8);
plt.title(title_list8, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,Aug_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,9);
plt.title(title_list9, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21  ## change this accordingly
colormesh_1 = spatial_figure(ax,Sep_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('SON',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,10);
plt.title(title_list10, fontsize = 30, y=1)
colormap=cmap;colorbar_min=0.5e20;colorbar_max=5e21 ## change this accordingly
colormesh_1 = spatial_figure(ax,Oct_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,11);
plt.title(title_list11, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21 ## change this accordingly
colormesh_1 = spatial_figure(ax,Nov_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

ax = plt.subplot(3,4,12);
plt.title(title_list12, fontsize = 30, y=1)
colormap=cmap; colorbar_min=0.5e20;colorbar_max=5e21 ## change this accordingly
colormesh_1 = spatial_figure(ax,Dec_mon_mean_model_ratio,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('DJF',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	

cbar_ax = fig.add_axes([0.10, 0.02, 0.80, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'GC Ratio (molucules/cm2/kg)',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.25, hspace=0.05);
#plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'mon_mean_model_ratio_without_mask.png',bbox_inches='tight')   ##########	
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'mon_mean_model_ratio_with_ocean_mask.png',bbox_inches='tight') 
plt.show()
