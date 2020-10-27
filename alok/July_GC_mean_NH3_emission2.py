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
from iris.coords import DimCoord
from iris.cube import Cube
import iris

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
		#print (month, 'iN @ Month')
		if day<=9:
			day = '0'+str(day);
		else:
			day =str(day);
		#print (day, 'iN @ Day')
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
		#print (day, 'iN @ Day')
		emissions_data_files = '/scratch/uptrop/em440/for_Alok/gc_ncdf/emissions/HEMCO_diagnostics.2016'+str(month)+str(day)+'0000.nc'
		return emissions_data_files
	
	#Satellites daily files 
	sat_data_files = create_filename_sat(month,day)
	ncf_sat = nc4.Dataset(sat_data_files,mode='r')
	lat_sat = ncf_sat.variables['LAT'][:]
	lon_sat = ncf_sat.variables['LON'][:]
	nh3_GC_column = ncf_sat.variables['IJ-AVG-S__NH3'][:]     		#NH3 tracer 'ppbv'
	airdensity_sat = ncf_sat.variables['TIME-SER__AIRDEN'][:]	#Air density 'molecules/cm3'
	bxheight_sat = ncf_sat.variables['BXHGHT-S__BXHEIGHT'][:]	#Grid Box height 'm'
	#print (nh3_GC_column.shape, airdensity_sat.shape, bxheight_sat.shape)
	
	bxheight1_sat = bxheight_sat*100 #Grid Box height 'cm'
	airdensity1_sat = airdensity_sat * bxheight1_sat # Air density 'molecules/cm3' * Grid Box height 'cm' = molecules/cm2
		
	nh3_GC_column_A = (nh3_GC_column/1e9)*airdensity1_sat#molucules/cm2
	#print (nh3_GC_column_A.shape)
	
	nh3_GC_column_B = np.nansum(nh3_GC_column_A, axis=0) #sum over all model vertical layers #molucules/cm2
	#print (nh3_GC_column_B.shape)
	
	#nh3_GC_column_C = nh3_GC_column_B/NA  #unit moles (NH3) /cm2
	#print (nh3_GC_column_C.shape)
	
	#emissions daily files
	emission_data_files = create_filename_emission(month,day)
	
	ncf_emission = nc4.Dataset(emission_data_files,mode='r')
	lat_emission = ncf_emission.variables['lat'][:]
	lon_emission = ncf_emission.variables['lon'][:]
	#area_emission = ncf_emission.variables['AREA'][:]  				#unit m2
	nh3_emission_total = ncf_emission.variables['EmisNH3_Total'][:]  	#Unit kg
	#nh3_emission_anthro = ncf_emission.variables['EmisNH3_Anthro'][:] 	#Unit kg
	#print (area_emission.shape, 'area_emission.shape')
	#print (nh3_emission_total.shape, 'nh3_emission_total.shape')
	
	nh3_emission_total_A = nh3_emission_total[0,:,:,:]
	#print (nh3_emission_total_A.shape, 'nh3_emission_total_A.shape')
	#sum over all vertical layers
	nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)			#kg  --> sum over all vertical layers
	#print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 		#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)           	#moles(NH3)/cm2
	#print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	#nh3_emission_total_E = nh3_emission_total_D*NA                      	#molecules/cm2
	#print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')
	#print (lat_GC_column.shape, lon_GC_column.shape, 'lat lon shape of sat file')
	#print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')
	return lat_sat,lon_sat,nh3_GC_column_B, lat_emission,lon_emission,nh3_emission_total_B

def monthly_mean_cal():	
	months=np.arange(1,13)
	time_series = np.empty(len(months));	

	GC_column_mon_mean = np.empty((len(time_series),115,177))
	GC_column_mon_mean[:] = np.nan
	
	emission_mon_mean = np.empty((len(time_series),115,177))
	emission_mon_mean[:] = np.nan
	
	for imonth in months:
		GC_column_nh3_mon_mean = np.empty((115,177))
		GC_column_nh3_mon_mean[:] = np.nan
		
		emission_nh3_mon_mean = np.empty((115,177))
		emission_nh3_mon_mean[:] = np.nan
				
		days = days_months(imonth)
		for iday in days:	
			lat_sat,lon_sat,nh3_GC_column_B, lat_emission,lon_emission,nh3_emission_total_B = daily_data(imonth,iday)
			GC_column_nh3_mon_mean = np.dstack((GC_column_nh3_mon_mean,nh3_GC_column_B))
			emission_nh3_mon_mean = np.dstack((emission_nh3_mon_mean,nh3_emission_total_B))
			
			#print (GC_column_nh3_mon_mean.shape)
		GC_column_mon_mean[imonth-1,:,:] = np.nanmean(GC_column_nh3_mon_mean,axis=2)
		#emission_mon_mean[imonth-1,:,:] = np.nanmean(emission_nh3_mon_mean,axis=2)
		emission_mon_mean[imonth-1,:,:] = np.nansum(emission_nh3_mon_mean,axis=2)
	#print (lat_sat.shape, lon_sat.shape, 'lat lon shape of sat file')
	#print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')	
	
	return time_series, lat_sat, lon_sat, emission_mon_mean, GC_column_mon_mean 

time_series, lat, lon, emission_mon_mean, GC_column_mon_mean = monthly_mean_cal()
#print(GC_column_mon_mean.shape,'!GC_column_mon_mean.shape')
#print(emission_mon_mean.shape,'!emission_mon_mean.shape')

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
#print (time)

cube1 = Cube(GC_column_mon_mean,
					dim_coords_and_dims=[(latitude, 1),
										(longitude, 2),
										(time, 0)])

cube2 = Cube(emission_mon_mean,
					dim_coords_and_dims=[(latitude, 1),
										(longitude, 2),
										(time, 0)])

regridded_data_GC_column = cube1.interpolate([('latitude', lat01), ('longitude', lon01)],
                           iris.analysis.Linear())
print(regridded_data_GC_column.shape,'regrid_GC_column_shape')

regridded_data_emission = cube2.interpolate([('latitude', lat01), ('longitude', lon01)],
                           iris.analysis.Linear())
print(regridded_data_emission.shape, 'regrid_emission_shape')



lat_n = regridded_data_emission.coord('latitude')
#print(lat_n.shape)
lat_n =lat_n.points[:]
lon_n = regridded_data_emission.coord('longitude')
#print(lon_n.shape)
lon_n =lon_n.points[:]

lat_n_min,lon_n_min = np.nanmin(lat_n),np.nanmin(lon_n)
lat_n_max,lon_n_max = np.nanmax(lat_n),np.nanmax(lon_n)
#print (lat_n_min, 'lat_min_regridded_data')
#print (lon_n_min, 'lon_min_regridded_data')
#print (lat_n_max, 'lat_max_regridded_data')
#print (lon_n_max, 'lon_max_regridded_data')

lat_gc_uk = lat_n[172:279]
#print (lat_gc_uk.shape, lat_gc_uk, 'lat_gc_uk_shape')
lon_gc_uk = lon_n[50:176]
#print (lon_gc_uk.shape, lon_gc_uk, 'lon_gc_uk_shape')



regridded_data_GC_column = regridded_data_GC_column[:].data
regridded_data_GC_column_uk = regridded_data_GC_column[:,172:279,50:176]
print (regridded_data_GC_column_uk.shape, 'regridded_data_GC_column_uk.shape')

July_GC_NH3_column = regridded_data_GC_column_uk[6, :, :]

title_listA = 'JULY GEOS-Chem NH3 Column'	

# Plot Figure 
fig = plt.figure(facecolor='White',figsize=[11,12.5]);pad= 1.1; 
#plt.suptitle(title_listC, fontsize = 35, y=1.001)

ax = plt.subplot(1,1,1);
plt.title(title_listA, fontsize =20, y=1)
colormap=cmap; colorbar_min=0;colorbar_max=3e15## change this accordingly
colormesh_1 = spatial_figure(ax,July_GC_NH3_column,lon_gc_uk,lat_gc_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
#ax.annotate('NAEI_derived = {0:.2f}'.format(NAEI_total_emission),xy=(0.3,0.001), xytext=(0, pad),
#		xycoords='axes fraction', textcoords='offset points',
#		ha='center', va='bottom',rotation='horizontal',fontsize=15)



cbar_ax = fig.add_axes([0.10, 0.02, 0.80, 0.01])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'JULY GC NH3 Column (molecules/cm2) ',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.25, hspace=0.05);
#plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'mon_mean_model_ratio_regrid_without_mask.png',bbox_inches='tight')   ##########	
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'July_GC_mean_NH3_column2.png',bbox_inches='tight')
plt.show()

