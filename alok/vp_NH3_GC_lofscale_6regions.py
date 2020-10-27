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
from netCDF4 import Dataset, num2date
from matplotlib.colors import LogNorm
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d
# cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')
cmap = cm.jet

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

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area

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

def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
	"""
	fill the region outside the box with 0
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	if (colum_s> colum_e):
		cache = colum_e; colum_e = colum_s; colum_s = cache;
	if (row_s> row_e):
		cache = row_e; row_e = row_s; row_s = cache;
	mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
	# plt.imshow(mask,origin='lower');plt.show()
	mask[0:row_s,:] =0; mask[row_e:-1,:] =0
	# plt.imshow(mask,origin='lower');plt.show()
	return mask	

def mask_weight(region_key,lon,lat,return_option,reverse=False):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	crop the sepecified region eithr as a box or an administrtative polygon
	input: 
		region_ky: region name, say, India
		lon and lat of your data
	output: depent on output_option
		if output_option == 'mask': output mask (1 for mask and nan for others)
		elif output_option == 'area': output area of a mask
		elif output_option == 'area_weight': output weight of area against the total area of the mask, this is useful when you do an area-weighted mean
	"""
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	print (lon_res, 'lon_res')
	lons,lats = np.meshgrid(lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
	print (area, 'area')

	##OCEAN_MASKS FOR COUNTRIES
	# ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	ocean_mask = sio.loadmat('/scratch/uptrop/ap744/shapefiles/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	## define your regions here
	#box_region_dic={'North_India':[75,80,25,30],'South_India':[75,80,15,20],'East_China':[108,123,30,40],'West_China':[85,100,30,40],'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	box_region_dic={'Wales':[-4,-3,51.5,53],'South_England':[-2.5,-1,51,52],'East_England':[0,1.5,51.5,53],'NE_England':[-2.5,-1.5,54.5,55.5],'N_Ireland':[-8,-6,54,55],'Scotland':[-5.5,-3,56.25,57.5],'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	
	if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
		mask= ocean_mask[region_key][:]
	elif  region_key in box_region_dic:
		mask= ocean_mask['All'][:]
		box = box_region_dic[region_key]
		mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
	else:
		print ("error region name")
	
	# interpolate from 360*720 to your grids
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
	# plt.imshow(mask,origin='lower');plt.show()
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	# weight each grid cell by its area weight against the total area
	if reverse:    ## note this Flase by default, but can be used to exclude the specified region from a larger map
		mask=1-mask
	mask[mask==0] = np.nan
	grid_area=np.multiply(mask,area); 
	mask_weighted = np.divide(grid_area,np.nansum(np.nansum(grid_area,axis=1),axis=0))
	if return_option == 'mask': return mask
	elif return_option == 'area': return grid_area
	elif return_option == 'area_weight': return mask_weighted

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
		#print (month, 'iN @ Month')
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
	nh3_sat = ncf_sat.variables['IJ-AVG-S__NH3'][:]     		#NH3 tracer 'ppbv'
	airdensity_sat = ncf_sat.variables['TIME-SER__AIRDEN'][:]	#Air density 'molecules/cm3'
	bxheight_sat = ncf_sat.variables['BXHGHT-S__BXHEIGHT'][:]	#Grid Box height 'm'
	#print (nh3_sat.shape, 'nh3_shape')
	
	#bxheight1_sat = bxheight_sat*100 #Grid Box height 'cm'
	#airdensity1_sat = airdensity_sat * bxheight1_sat # Air density 'molecules/cm3' * Grid Box height 'cm' = molecules/cm2
	#nh3_sat_A = nh3_sat*airdensity1_sat #molucules/cm2
	#print (nh3_sat_A.shape)
	#nh3_sat_B = np.nansum(nh3_sat_A, axis=0) #sum over all model vertical layers
	#print (nh3_sat_B.shape)
	
	#emissions daily files
	emission_data_files = create_filename_emission(month,day)
	
	ncf_emission = nc4.Dataset(emission_data_files,mode='r')
	lat_emission = ncf_emission.variables['lat'][:]
	lon_emission = ncf_emission.variables['lon'][:]
	area_emission = ncf_emission.variables['AREA'][:]  					#unit m2
	nh3_emission_total = ncf_emission.variables['EmisNH3_Total'][:]  	#Unit kg
	nh3_emission_anthro = ncf_emission.variables['EmisNH3_Anthro'][:] 	#Unit kg
	#print (area_emission.shape, 'area_emission.shape')
	#print (nh3_emission_total.shape, 'nh3_emission_total.shape')
	
	nh3_emission_total_A = nh3_emission_total[0,:,:,:]
	#print (nh3_emission_total_A.shape, 'nh3_emission_total_A.shape')
	#sum over all vertical layers
	#nh3_emission_total_B = np.nansum(nh3_emission_total_A, axis=0)			#kg  --> sum over all vertical layers
	#print (nh3_emission_total_B.shape, 'nh3_emission_total_B.shape')
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 		#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)           	#moles(NO2)/cm2
	#print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	#nh3_emission_total_E = nh3_emission_total_D*NA                      	#molecules/cm2
	#print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')
	#print (lat_sat.shape, lon_sat.shape, 'lat lon shape of sat file')
	#print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')
	return lat_sat,lon_sat,nh3_sat, lat_emission,lon_emission,nh3_emission_total_A

def monthly_mean_cal():	
	months=np.arange(1,3)
	time_series = np.empty(len(months));	

	sat_mon_mean = np.empty((len(time_series),40,115,177))
	sat_mon_mean[:] = np.nan
	
	emission_mon_mean = np.empty((len(time_series),47,115,177))
	emission_mon_mean[:] = np.nan
	
	for imonth in months:
		sat_nh3_mon_mean = np.empty((40,115,177))
		sat_nh3_mon_mean[:] = np.nan
		
		emission_nh3_mon_mean = np.empty((47,115,177))
		emission_nh3_mon_mean[:] = np.nan
				
		days = days_months(imonth)
		new_sat_nh3_mon_mean = np.empty(((len(days)+1),40,115,177))
		new_sat_nh3_mon_mean[:] = np.nan
		
		new_emission_nh3_mon_mean = np.empty(((len(days)+1),47,115,177))
		new_emission_nh3_mon_mean[:] = np.nan
		
		for iday in days:	
			#print (iday, '!!!iday!!!')
			lat_sat,lon_sat,nh3_sat, lat_emission,lon_emission,nh3_emission_total_A = daily_data(imonth,iday)
			#print (nh3_sat.shape, 'nh3_sat')
			sat_nh3_mon_mean = np.concatenate((sat_nh3_mon_mean,nh3_sat), axis=0)
			#print (sat_nh3_mon_mean.shape , 'sat_nh3_mon_mean.shape -concat')
			new_sat_nh3_mon_mean[iday] = sat_nh3_mon_mean[iday*40:iday*40+40]
			# plt.plot(new_sat_nh3_mon_mean[i,:])
			#print(new_sat_nh3_mon_mean.shape)
			
			emission_nh3_mon_mean = np.concatenate((emission_nh3_mon_mean,nh3_emission_total_A), axis=0)
			new_emission_nh3_mon_mean[iday] = emission_nh3_mon_mean[iday*47:iday*47+47]
			
		sat_mon_mean[imonth-1,:,:,:] = np.nanmean(new_sat_nh3_mon_mean,axis=0)
		emission_mon_mean[imonth-1,:,:,:] = np.nanmean(new_emission_nh3_mon_mean,axis=0)
	#print (lat_sat.shape, lon_sat.shape, 'lat lon shape of sat file')
	#print (lat_emission.shape, lon_emission.shape, 'lat lon shape of emission file')	
	
	return time_series, lat_sat, lon_sat, emission_mon_mean, sat_mon_mean 

time_series, lat, lon, emission_mon_mean, sat_mon_mean = monthly_mean_cal()
#print(sat_mon_mean.shape,'!sat_mon_mean.shape')
#print(emission_mon_mean.shape,'!emission_mon_mean.shape')

def Wales_ts_mon():
	region_key='Wales';return_option='area_weight'
	mask_walesA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	wales_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_walesA,sat_mon_mean),axis=3),axis=2)
	#print (wales_NH3_GC.shape)
	
	return wales_NH3_GC

def NE_England_ts_mon():
	region_key='NE_England';return_option='area_weight'
	mask_NE_EnglandA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	NE_England_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_NE_EnglandA,sat_mon_mean),axis=3),axis=2)
	#print (NE_England_NH3_GC.shape)
	
	return NE_England_NH3_GC

def South_England_ts_mon():
	region_key='South_England';return_option='area_weight'
	mask_South_EnglandA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	South_England_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_South_EnglandA,sat_mon_mean),axis=3),axis=2)
	#print (South_England_NH3_GC.shape)
	
	return South_England_NH3_GC

def East_England_ts_mon():
	region_key='East_England';return_option='area_weight'
	mask_East_EnglandA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	East_England_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_East_EnglandA,sat_mon_mean),axis=3),axis=2)
	print (East_England_NH3_GC)
	
	return East_England_NH3_GC

def N_Ireland_ts_mon():
	region_key='N_Ireland';return_option='area_weight'
	mask_N_IrelandA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	N_Ireland_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_N_IrelandA,sat_mon_mean),axis=3),axis=2)
	print (N_Ireland_NH3_GC)
	
	return N_Ireland_NH3_GC

def Scotland_ts_mon():
	region_key='Scotland';return_option='area_weight'
	mask_ScotlandA=mask_weight(region_key,lon.copy(),lat.copy(),return_option,reverse=False);
	## weighted average of your data
	Scotland_NH3_GC 		= np.nansum(np.nansum(np.multiply(mask_ScotlandA,sat_mon_mean),axis=3),axis=2)
	#print (Scotland_NH3_GC.shape)
	
	return Scotland_NH3_GC

wales_NH3_GC = Wales_ts_mon()
#print (wales_NH3_GC.shape, 'wales_NH3_GC.shape' )
wales_NH3_GC_vp1 = wales_NH3_GC[:,:30]
wales_NH3_GC_vp =wales_NH3_GC_vp1.transpose()
#print (wales_NH3_GC_vp.shape, 'wales_NH3_GC_vp.shape')

NE_England_NH3_GC = NE_England_ts_mon()
NE_England_NH3_GC_vp1 = NE_England_NH3_GC[:,:30]
NE_England_NH3_GC_vp =NE_England_NH3_GC_vp1.transpose()

South_England_NH3_GC = South_England_ts_mon()
South_England_NH3_GC_vp1 = South_England_NH3_GC[:,:30]
South_England_NH3_GC_vp =South_England_NH3_GC_vp1.transpose()

East_England_NH3_GC = East_England_ts_mon()
East_England_NH3_GC_vp1 = East_England_NH3_GC[:,:30]
East_England_NH3_GC_vp =East_England_NH3_GC_vp1.transpose()

N_Ireland_NH3_GC = N_Ireland_ts_mon()
N_Ireland_NH3_GC_vp1 = N_Ireland_NH3_GC[:,:30]
N_Ireland_NH3_GC_vp =N_Ireland_NH3_GC_vp1.transpose()

Scotland_NH3_GC = Scotland_ts_mon()
Scotland_NH3_GC_vp1 = Scotland_NH3_GC[:,:30]
Scotland_NH3_GC_vp =Scotland_NH3_GC_vp1.transpose()

"""
Y = [58, 189, 320, 454, 589, 726, 864, 1004, 1146, 1290, 
	 1436, 1584, 1759, 1988, 2249, 2517, 2792, 3074, 3439, 3896,
     4375, 4879, 5413, 5980, 6585, 7237, 7943, 8846, 9936, 11021]
     #6820, 7253, 7700, 8160, 8633, 9120, 9620, 10133, 10660, 11200,
     #11754, 12321, 12901, 13495, 14102, 14724, 15359, 16009, 16673, 17352],
     #18046, 18757, 19484, 20229, 20993, 21777, 22582, 23412, 24268, 25153,
     #26071, 27024, 28018, 29058, 30150, 31301, 32518, 33811, 35190, 36666,
     #38254	,39968, 41825, 43844, 46046, 48456, 51099, 54006, 57210, 60747]# 64657, 68986, 73782, 79100,	85000]
# print (Y)

X= np.linspace(1,12,12)
#print (X)

#levels = [0, 1.75, 3.50, 5.25, 7.0, 8.75, 10.5, 12.25, 14.0, 15.75, 17.50]
#levels = np.linspace(0,19,11)
levels = np.logspace(-2.5, 0.6, 11)
#print (levels)
# tick_levels = [0.001, 0.002, 0.007, 0.02, 0.05, 0.14, 0.38, 1.02, 2.75, 7.41, 19.95]
tick_levels =[0.003, 0.006, 0.013, 0.027, 0.055, 0.112, 0.229,0.467,0.954,1.949,3.981]

fig = plt.figure(facecolor='White',figsize=[16,11]);pad= 0.5; 
# plt.suptitle ('Vertical Profile NO$_2$ ', fontsize = 32, y=0.95)

ax = fig.add_subplot(321)
plt.contourf(wales_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(wales_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(wales_NH3_GC_vp, cmap = cm.jet)
plt.title('Wales',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9,11]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16','Dec-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
ax.axes.get_xaxis().set_visible(False)

ax = fig.add_subplot(322)
plt.contourf(NE_England_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(wales_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(wales_NH3_GC_vp, cmap = cm.jet)
plt.title('NE_England',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9,11]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16','Dec-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)

ax = fig.add_subplot(323)
plt.contourf(South_England_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(South_England_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(South_England_NH3_GC_vp, cmap = cm.jet)
plt.title('South_England',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
#ax.axes.get_xaxis().set_visible(False)

ax = fig.add_subplot(324)
plt.contourf(East_England_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(East_England_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(wales_NH3_GC_vp, cmap = cm.jet)
plt.title('East_England',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
ax.axes.get_yaxis().set_visible(False)

ax = fig.add_subplot(325)
plt.contourf(N_Ireland_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(East_England_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(wales_NH3_GC_vp, cmap = cm.jet)
plt.title('N_Ireland',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
ax.axes.get_yaxis().set_visible(False)

ax = fig.add_subplot(326)
plt.contourf(Scotland_NH3_GC_vp, levels, norm=LogNorm(), cmap = cm.jet,extend='both')
#plt.contourf(East_England_NH3_GC_vp, levels, cmap = cm.jet)
#plt.contourf(wales_NH3_GC_vp, cmap = cm.jet)
plt.title('Scotland',fontsize = 30, y=1.01)
plt.ylabel('Height (meters)', fontsize = 30, y =0.5)
plt.xlabel('Months', fontsize = 25, y=0.01)
cb1 = plt.colorbar(aspect = 25, ticks=levels, orientation='vertical',extend='both',fraction=0.1,pad=0.01)
cb1.set_label('NH$_3$  (ppbv)', fontsize = 25)
cb1.ax.set_yticklabels(tick_levels)
cb1.ax.tick_params(labelsize=25)
tick_locs = [0,3,6,9]
tick_lbls = ['Jan-16','Apr-16', 'Jul-16','Oct-16']
plt.xticks(tick_locs, tick_lbls, fontsize = 25)
tick_locs1 = [1,6,12,18,24,30]
tick_lbls1 = ['58','726', '1584','3074', '5980','11021']
plt.yticks(tick_locs1, tick_lbls1, fontsize = 30, y=1.01)
ax.axes.get_yaxis().set_visible(False)

plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.21, hspace=0.11);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'vertical_profile_logscale_6_boxes.png',bbox_inches='tight')   ##########	
plt.show()"""
