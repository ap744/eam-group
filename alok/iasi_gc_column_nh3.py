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
	
	nh3_sat_A = (nh3_sat*airdensity1_sat)/1e9 #molucules/cm2
	print (nh3_sat_A.shape)
	
	nh3_sat_B = np.nansum(nh3_sat_A, axis=0) #sum over all model vertical layers #molucules/cm2
	print (nh3_sat_B.shape)
	
	#nh3_sat_C = nh3_sat_B/NA  #unit moles (NH3) /cm2
	#print (nh3_sat_C.shape)
	
	#emissions daily files
	emission_data_files = create_filename_emission(month,day)
	
	ncf_emission = nc4.Dataset(emission_data_files,mode='r')
	lat_emission = ncf_emission.variables['lat'][:]
	lon_emission = ncf_emission.variables['lon'][:]
	#area_emission = ncf_emission.variables['AREA'][:]  					#unit m2
	nh3_emission_total = ncf_emission.variables['EmisNH3_Total'][:]  	#Unit kg
	#nh3_emission_anthro = ncf_emission.variables['EmisNH3_Anthro'][:] 	#Unit kg
	#print (area_emission.shape, 'area_emission.shape')
	print (nh3_emission_total.shape, 'nh3_emission_total.shape')
	
	nh3_emission_total_A = nh3_emission_total[0,:,:,:]
	print (nh3_emission_total_A.shape, 'nh3_emission_total_A.shape')
	#sum over all vertical layers
	nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)			#kg  --> sum over all vertical layers
	print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 		#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)           	#moles(NH3)/cm2
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

cube1 = Cube(sat_mon_mean,
					dim_coords_and_dims=[(latitude, 1),
										(longitude, 2),
										(time, 0)])

cube2 = Cube(emission_mon_mean,
					dim_coords_and_dims=[(latitude, 1),
										(longitude, 2),
										(time, 0)])

regridded_data_sat = cube1.interpolate([('latitude', lat01), ('longitude', lon01)],
                           iris.analysis.Linear())
print(regridded_data_sat.shape,'regrid_sat')

regridded_data_emission = cube2.interpolate([('latitude', lat01), ('longitude', lon01)],
                           iris.analysis.Linear())
print(regridded_data_emission.shape, 'regrid_emission')



lat_n = regridded_data_emission.coord('latitude')
#print(lat_n.shape)
lat_n =lat_n.points[:]
lon_n = regridded_data_emission.coord('longitude')
#print(lon_n.shape)
lon_n =lon_n.points[:]

lat_n_min,lon_n_min = np.nanmin(lat_n),np.nanmin(lon_n)
lat_n_max,lon_n_max = np.nanmax(lat_n),np.nanmax(lon_n)
print (lat_n_min, 'lat_min_regridded_data')
print (lon_n_min, 'lon_min_regridded_data')
print (lat_n_max, 'lat_max_regridded_data')
print (lon_n_max, 'lon_max_regridded_data')

lat_gc_uk = lat_n[172:279]
print (lat_gc_uk.shape, lat_gc_uk, 'lat_gc_uk_shape')
lon_gc_uk = lon_n[50:176]
print (lon_gc_uk.shape, lon_gc_uk, 'lon_gc_uk_shape')

regridded_data_emission = regridded_data_emission[:].data
regridded_data_emission_uk_iasi = regridded_data_emission[:,172:279,50:176].copy()
regridded_data_emission_uk_naei = regridded_data_emission[:,172:279,50:176].copy()
print (regridded_data_emission_uk_iasi.shape , 'regridded_data_emission_uk_iasi.shape')

#mask low emission mainly over ocean
regridded_data_emission_uk_iasi[regridded_data_emission_uk_iasi <= 5] = np.nan

regridded_data_sat = regridded_data_sat[:].data
regridded_data_sat_uk = regridded_data_sat[:,172:279,50:176]
print (regridded_data_sat_uk.shape, 'regridded_data_sat_uk.shape')

regridded_data_sat_uk_annual = np.nansum(regridded_data_sat_uk,axis=0)/1e15
print (regridded_data_sat_uk_annual,'regridded_data_sat_uk_annual')

regridded_data_sat_uk_annual_total = np.nansum(np.nansum(regridded_data_sat_uk_annual,axis=1),axis=0)
print (regridded_data_sat_uk_annual_total,'regridded_data_sat_uk_annual_total')


#calculation of GC ratio (NH3 column concentration / NH3 emission)
Jan_mon_mean_model_ratio = regridded_data_sat_uk[0,:,:]/regridded_data_emission_uk_iasi[0,:,:]
Feb_mon_mean_model_ratio = regridded_data_sat_uk[1,:,:]/regridded_data_emission_uk_iasi[1,:,:]
Mar_mon_mean_model_ratio = regridded_data_sat_uk[2,:,:]/regridded_data_emission_uk_iasi[2,:,:]
Apr_mon_mean_model_ratio = regridded_data_sat_uk[3,:,:]/regridded_data_emission_uk_iasi[3,:,:]
May_mon_mean_model_ratio = regridded_data_sat_uk[4,:,:]/regridded_data_emission_uk_iasi[4,:,:]
Jun_mon_mean_model_ratio = regridded_data_sat_uk[5,:,:]/regridded_data_emission_uk_iasi[5,:,:]
Jul_mon_mean_model_ratio = regridded_data_sat_uk[6,:,:]/regridded_data_emission_uk_iasi[6,:,:]
Aug_mon_mean_model_ratio = regridded_data_sat_uk[7,:,:]/regridded_data_emission_uk_iasi[7,:,:]
Sep_mon_mean_model_ratio = regridded_data_sat_uk[8,:,:]/regridded_data_emission_uk_iasi[8,:,:]
Oct_mon_mean_model_ratio = regridded_data_sat_uk[9,:,:]/regridded_data_emission_uk_iasi[9,:,:]
Nov_mon_mean_model_ratio = regridded_data_sat_uk[10,:,:]/regridded_data_emission_uk_iasi[10,:,:]
Dec_mon_mean_model_ratio = regridded_data_sat_uk[11,:,:]/regridded_data_emission_uk_iasi[11,:,:]
print (Jan_mon_mean_model_ratio.shape,'Jan_model_ratio.shape')

#Reading IASI column concentration
iasi_nh3_file = nc4.Dataset('/scratch/uptrop/em440/for_Alok/iasi_ncdf/iasi_nh3_uk_oversampled_2008-2018_0.1_jul2020.nc',mode='r')
lat_iasi = iasi_nh3_file.variables['lat'][:]
lon_iasi = iasi_nh3_file.variables['lon'][:]
iasi_nh3 = iasi_nh3_file.variables['iasi_nh3'][:] #unit molecules/cm2
lat_iasi_min,lon_iasi_min = np.nanmin(lat_iasi),np.nanmin(lon_iasi)
lat_iasi_max,lon_iasi_max = np.nanmax(lat_iasi),np.nanmax(lon_iasi)
print (lat_iasi_min, 'lat_min_iasi')
print (lon_iasi_min, 'lon_min_iasi')
print (lat_iasi_max, 'lat_max_iasi')
print (lon_iasi_max, 'lon_max_iasi')

lat_iasi_uk = lat_iasi[1:108]
print (lat_iasi_uk.shape, lat_iasi_uk, 'lat_iasi_uk_shape')
lon_iasi_uk = lon_iasi[1:127]
print (lon_iasi_uk.shape, lon_iasi_uk, 'lon_iasi_uk_shape')

print (iasi_nh3.shape, 'iasi_nh3.shape')
iasi_nh3_uk = iasi_nh3[:,1:108,1:127]
print (iasi_nh3_uk.shape, 'iasi_nh3_uk.shape')

iasi_nh3_uk_annual = np.nansum(iasi_nh3_uk,axis=0)/1e15
print (iasi_nh3_uk_annual,'iasi_nh3_uk_annual')

iasi_nh3_uk_annual_total = np.nansum(np.nansum(iasi_nh3_uk_annual,axis=1),axis=0)
print (iasi_nh3_uk_annual_total,'iasi_nh3_uk_annual_total')

################################# IASI MONTHLY EMISSION ##############################################
Jan_iasi_nh3_emission =  iasi_nh3_uk[0,:,:] / Jan_mon_mean_model_ratio
Feb_iasi_nh3_emission =  iasi_nh3_uk[1,:,:] / Feb_mon_mean_model_ratio
Mar_iasi_nh3_emission =  iasi_nh3_uk[2,:,:] / Mar_mon_mean_model_ratio
Apr_iasi_nh3_emission =  iasi_nh3_uk[3,:,:] / Apr_mon_mean_model_ratio
May_iasi_nh3_emission =  iasi_nh3_uk[4,:,:] / May_mon_mean_model_ratio
Jun_iasi_nh3_emission =  iasi_nh3_uk[5,:,:] / Jun_mon_mean_model_ratio
Jul_iasi_nh3_emission =  iasi_nh3_uk[6,:,:] / Jul_mon_mean_model_ratio
Aug_iasi_nh3_emission =  iasi_nh3_uk[7,:,:] / Aug_mon_mean_model_ratio 
Sep_iasi_nh3_emission =  iasi_nh3_uk[8,:,:] / Sep_mon_mean_model_ratio 
Oct_iasi_nh3_emission =  iasi_nh3_uk[9,:,:] / Oct_mon_mean_model_ratio
Nov_iasi_nh3_emission =  iasi_nh3_uk[10,:,:] / Nov_mon_mean_model_ratio
Dec_iasi_nh3_emission =  iasi_nh3_uk[11,:,:] / Dec_mon_mean_model_ratio
print (Jan_iasi_nh3_emission.shape, 'Jan_iasi_nh3_emission.shape')
#print (np.nanmax(Apr_iasi_nh3_emission), 'Apr_Max', np.nanmin(Apr_iasi_nh3_emission), 'Apr_Min')



#################################################################################################
################################################ NAEI EMISSION ##################################
#################################################################################################

#calculating annual emission
annual_emission = np.nansum(regridded_data_emission_uk_naei.data[:], axis=0)

#monthly scale factor
Jan_scale_factor = regridded_data_emission_uk_naei[0,:,:]/annual_emission #unit - unitless(kg/kg from GC model)
Feb_scale_factor = regridded_data_emission_uk_naei[1,:,:]/annual_emission
Mar_scale_factor = regridded_data_emission_uk_naei[2,:,:]/annual_emission
Apr_scale_factor = regridded_data_emission_uk_naei[3,:,:]/annual_emission
May_scale_factor = regridded_data_emission_uk_naei[4,:,:]/annual_emission
Jun_scale_factor = regridded_data_emission_uk_naei[5,:,:]/annual_emission
Jul_scale_factor = regridded_data_emission_uk_naei[6,:,:]/annual_emission
Aug_scale_factor = regridded_data_emission_uk_naei[7,:,:]/annual_emission
Sep_scale_factor = regridded_data_emission_uk_naei[8,:,:]/annual_emission
Oct_scale_factor = regridded_data_emission_uk_naei[9,:,:]/annual_emission
Nov_scale_factor = regridded_data_emission_uk_naei[10,:,:]/annual_emission
Dec_scale_factor = regridded_data_emission_uk_naei[11,:,:]/annual_emission
print (Jan_scale_factor.shape, 'Jan_scale_factor.shape')


#Reading NAEI emission data
naei_nh3_file = nc4.Dataset('/scratch/uptrop/em440/for_Alok/naei_nh3/NAEI_total_NH3_0.1x0.1_2016.nc',mode='r')
lat_naei = naei_nh3_file.variables['lat'][:]
lon_naei = naei_nh3_file.variables['lon'][:]
naei_nh3 = naei_nh3_file.variables['NH3'][:] 	#unit g/m2/yr
naei_area = naei_nh3_file.variables['area'][:] 	#unit m2
#naei_nh3[naei_nh3<0.006] = np.nan 


naei_nh3_area = (naei_nh3 * naei_area )/1000 # g/m2/yr * m2 = g/yr --> g/yr/1000 --->kg/yr
#naei_nh3_area_mon = naei_nh3_area/12 # kg/month
naei_nh3_area_mon = naei_nh3_area
naei_nh3_area_mon[naei_nh3_area_mon<50] =np.nan

lat_naei_min,lon_naei_min = np.nanmin(lat_naei),np.nanmin(lon_naei)
lat_naei_max,lon_naei_max = np.nanmax(lat_naei),np.nanmax(lon_naei)
print (lat_naei_min, 'lat_min_naei')
print (lon_naei_min, 'lon_min_naei')
print (lat_naei_max, 'lat_max_naei')
print (lon_naei_max, 'lon_max_naei')

lat_naei_uk = lat_naei[7:114]
print (lat_naei_uk.shape, lat_naei_uk, 'lat_naei_uk_shape')
lon_naei_uk = lon_naei[4:130]
print (lon_naei_uk.shape, lon_naei_uk, 'lon_naei_uk_shape')

print (naei_nh3_area_mon.shape, 'naei_nh3_area_mon.shape')
naei_nh3_uk = naei_nh3_area_mon[7:114,4:130]
print (naei_nh3_uk.shape, 'naei_nh3_uk.shape')

mask_UK = naei_nh3_uk.copy()
mask_UK[mask_UK>49] = 1

NAEI_total_emission_raw = np.nansum(np.nansum(naei_nh3_uk,axis=1),axis=0)
print (NAEI_total_emission_raw,'NAEI_total_emission_raw')

############################ NAEI MONTHLY EMISSION ######################################
Jan_naei_nh3_emission = Jan_scale_factor *  naei_nh3_uk
Feb_naei_nh3_emission = Feb_scale_factor *  naei_nh3_uk
Mar_naei_nh3_emission = Mar_scale_factor *  naei_nh3_uk
Apr_naei_nh3_emission = Apr_scale_factor *  naei_nh3_uk
May_naei_nh3_emission = May_scale_factor *  naei_nh3_uk
Jun_naei_nh3_emission = Jun_scale_factor *  naei_nh3_uk
Jul_naei_nh3_emission = Jul_scale_factor *  naei_nh3_uk
Aug_naei_nh3_emission = Aug_scale_factor *  naei_nh3_uk
Sep_naei_nh3_emission = Sep_scale_factor *  naei_nh3_uk
Oct_naei_nh3_emission = Oct_scale_factor *  naei_nh3_uk
Nov_naei_nh3_emission = Nov_scale_factor *  naei_nh3_uk
Dec_naei_nh3_emission = Dec_scale_factor *  naei_nh3_uk

print (Jan_naei_nh3_emission.shape, 'Jan_naei_nh3_emission.shape')

############################### Diff NAEI - IASI ########################################
diff_Jan_NAEI_IASI = Jan_naei_nh3_emission - Jan_iasi_nh3_emission
diff_Feb_NAEI_IASI = Feb_naei_nh3_emission - Feb_iasi_nh3_emission
diff_Mar_NAEI_IASI = Mar_naei_nh3_emission - Mar_iasi_nh3_emission
diff_Apr_NAEI_IASI = Apr_naei_nh3_emission - Apr_iasi_nh3_emission
diff_May_NAEI_IASI = May_naei_nh3_emission - May_iasi_nh3_emission
diff_Jun_NAEI_IASI = Jun_naei_nh3_emission - Jun_iasi_nh3_emission
diff_Jul_NAEI_IASI = Jul_naei_nh3_emission - Jul_iasi_nh3_emission
diff_Aug_NAEI_IASI = Aug_naei_nh3_emission - Aug_iasi_nh3_emission
diff_Sep_NAEI_IASI = Sep_naei_nh3_emission - Sep_iasi_nh3_emission
diff_Oct_NAEI_IASI = Oct_naei_nh3_emission - Oct_iasi_nh3_emission
diff_Nov_NAEI_IASI = Nov_naei_nh3_emission - Nov_iasi_nh3_emission
diff_Dec_NAEI_IASI = Dec_naei_nh3_emission - Dec_iasi_nh3_emission


Jan_Dec_NAEI = np.stack([Jan_naei_nh3_emission, Feb_naei_nh3_emission, Mar_naei_nh3_emission,Apr_naei_nh3_emission,May_naei_nh3_emission,Jun_naei_nh3_emission,Jul_naei_nh3_emission,Aug_naei_nh3_emission,Sep_naei_nh3_emission,Oct_naei_nh3_emission, Nov_naei_nh3_emission, Dec_naei_nh3_emission])
print (Jan_Dec_NAEI.shape)
Jan_Dec_NAEI_sum = np.nansum(Jan_Dec_NAEI, axis=0)
print (Jan_Dec_NAEI_sum.shape, 'naei sum')
Jan_Dec_NAEI_sum[Jan_Dec_NAEI_sum<50]=np.nan

NAEI_total_emission = np.nansum(np.nansum(Jan_Dec_NAEI_sum,axis=1),axis=0)
print (NAEI_total_emission,'NAEI_total_emission')

Jan_Dec_IASI = np.array([Jan_iasi_nh3_emission,Feb_iasi_nh3_emission,Mar_iasi_nh3_emission,Apr_iasi_nh3_emission,May_iasi_nh3_emission,Jun_iasi_nh3_emission,Jul_iasi_nh3_emission,Aug_iasi_nh3_emission,Sep_iasi_nh3_emission,Oct_iasi_nh3_emission,Nov_iasi_nh3_emission,Dec_iasi_nh3_emission])
print (Jan_Dec_IASI.shape)
Jan_Dec_IASI_sum = np.nansum(Jan_Dec_IASI, axis=0)
print (Jan_Dec_IASI_sum.shape, 'IASI sum')
#Jan_Dec_IASI_sum[Jan_Dec_IASI_sum<500]=np.nan

Jan_Dec_IASI_sum = Jan_Dec_IASI_sum*mask_UK
IASI_total_emission = np.nansum(np.nansum(Jan_Dec_IASI_sum,axis=1),axis=0)
print (IASI_total_emission,'IASI_total_emission')


title_listA = 'GC Column NH$_3$ (10$^{15}$molecules/cm$^2$)'	
title_listB = 'IASI Column NH$_3$ (10$^{15}$molecules/cm$^2$)'	


# Plot Figure 
fig = plt.figure(facecolor='White',figsize=[16,11]);pad= 1.1; 
#plt.suptitle(title_listC, fontsize = 35, y=1.001)

ax = plt.subplot(1,2,1);
plt.title(title_listA, fontsize = 25, y=1.03)
colormap=cmap; colorbar_min=1;colorbar_max=51## change this accordingly
colormesh_1 = spatial_figure(ax,regridded_data_sat_uk_annual,lon_naei_uk,lat_naei_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
#ax.annotate('GC_Column = {0:.2f}'.format(regridded_data_sat_uk_annual_total),xy=(0.4,0.001), xytext=(0, pad),
		#xycoords='axes fraction', textcoords='offset points',
		#ha='center', va='bottom',rotation='horizontal',fontsize=15)

ax = plt.subplot(1,2,2);
plt.title(title_listB, fontsize = 25, y=1.03)
colormap=cmap; colorbar_min=1;colorbar_max=51## change this accordingly
colormesh_1 = spatial_figure(ax,iasi_nh3_uk_annual,lon_naei_uk,lat_naei_uk,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False)
# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
		# xycoords='axes fraction', textcoords='offset points',
		# ha='center', va='bottom',rotation='horizontal',fontsize=30)	
#ax.annotate('IASI_Column = {0:.2f}'.format(iasi_nh3_uk_annual_total),xy=(0.4,0.001), xytext=(0, pad),
		#xycoords='axes fraction', textcoords='offset points',
		#ha='center', va='bottom',rotation='horizontal',fontsize=15)

cbar_ax = fig.add_axes([0.10, 0.08, 0.80, 0.025])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=20) 
cbar_ax.annotate(r'Annual NH$_3$ (10$^{15}$molecules/cm$^2$)',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.25, hspace=0.05);
#plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'mon_mean_model_ratio_regrid_without_mask.png',bbox_inches='tight')   ##########	
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'IASI_GC_Column_NH3.png',bbox_inches='tight')
plt.show()

