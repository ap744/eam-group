import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset, num2date
import datetime
import netCDF4 as nc4
import matplotlib.cm as cm
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.patches import Polygon
from matplotlib import gridspec
# cmap = cm.get_cmap('brewer_RdBu_11')
cmap = cm.jet
cmap = cm.rainbow
cmap = cm.YlOrRd
cmap1 = cm.coolwarm

Today_date=datetime.datetime.now().strftime("%Y%m%d")
title_list2  = 'NO2 Emissions trends'

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
	
def spatial_figure1(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=False): #c_bad,c_under,c_over,c_number=20,
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
	
	def draw_screen_poly( lats_XX, lons_XX, map):
		x, y = map( lons_XX, lats_XX)
		xy = zip(x,y)
		poly = Polygon( list(xy), facecolor='None', alpha=0.9, edgecolor='black',linewidth=5)
		plt.gca().add_patch(poly)
	lats_wales = [ 51.5, 53, 53, 51.5 ]
	lons_wales = [ -4, -4, -3, -3 ]

	lats_south_england = [ 51, 52, 52, 51 ]
	lons_south_england = [ -2.5, -2.5, -1, -1 ]

	lats_east_england = [ 51.5, 53, 53, 51.5 ]
	lons_east_england = [ 0, 0, 1.5, 1.5 ]

	lats_NE_england = [ 54.5, 55.5, 55.5, 54.5 ]
	lons_NE_england = [ -2.5, -2.5, -1.5, -1.5 ]
	
	lats_N_Ireland = [ 54, 55, 55, 54 ]
	lons_N_Ireland = [ -8, -8, -6, -6 ]
	
	lats_scotland = [ 56.25, 57.5, 57.5, 56.25 ]
	lons_scotland = [ -5.5, -5.5, -3, -3 ]
	
	

	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=30)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=30)
	# Add Coastlines, States, and Country Boundaries
	# map.drawcoastlines(); map.drawcountries() #map.drawstates(); # draw border lines, here only coast lines
	map.readshapefile('/scratch/uptrop/ap744/shapefiles/Shapfiles_india/World_shp/World','World', linewidth=2) # to add a shapefile on a map
	# map.readshapefile('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Shapfiles_india/World_shp/Export_Output','Export_Output')
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	draw_screen_poly( lats_wales, lons_wales, map )
	draw_screen_poly( lats_south_england, lons_south_england, map )
	draw_screen_poly( lats_east_england, lons_east_england, map )
	draw_screen_poly( lats_NE_england, lons_NE_england, map )
	draw_screen_poly( lats_N_Ireland, lons_N_Ireland, map )
	draw_screen_poly( lats_scotland, lons_scotland, map )

	#masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)   #  use 20 color bins, this can be changed
	cmap.set_bad([1,1,1],alpha = 1.0);
	if bad_data:
		cmap.set_under('w');cmap.set_over('k')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	return colormesh	
	
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
	lons,lats = np.meshgrid(lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

	##OCEAN_MASKS FOR COUNTRIES
	# ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	ocean_mask = sio.loadmat('/exports/csce/datastore/geos/groups/alok_pandey/Research_Work/Mat_File/Euro_USA_AUS_BRICS_STA_720_360.mat')  ## change this accordingly
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	## define your regions here
	box_region_dic={'North_India':[75,80,25,30],'South_India':[75,80,15,20],'East_China':[110,125,30,40],'West_China':[85,100,30,40],'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
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
	if reverse:    ## note this Flase by default, but can be used to exclude the specified  region from a larger map
		mask=1-mask
	mask[mask==0] = np.nan
	grid_area=np.multiply(mask,area); 
	mask_weighted = np.divide(grid_area,np.nansum(np.nansum(grid_area,axis=1),axis=0))
	if return_option == 'mask': return mask
	elif return_option == 'area': return grid_area
	elif return_option == 'area_weight': return mask_weighted
	

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
	
	nh3_emission_total_C = nh3_emission_total_B/(area_emission) 			#kg(NH3)/m2
	nh3_emission_total_E = nh3_emission_total_C*1000 						#g/m2/month as calculating fora a month
	
	##below 6 lines to convert kg to molecules/cm2
	#nh3_emission_total_C = nh3_emission_total_B/(area_emission*1.0e4) 	#kg(NH3)/cm2
	#print (nh3_emission_total_C.shape, 'nh3_emission_total_C.shape')
	#nh3_emission_total_D = nh3_emission_total_C/(mNH3*1.0e-3)          #moles(NO2)/cm2
	#print (nh3_emission_total_D.shape, 'nh3_emission_total_D.shape')
	#nh3_emission_total_E = nh3_emission_total_D*NA                     #molecules/cm2
	#print (nh3_emission_total_E.shape, 'nh3_emission_total_E.shape')	
	
	return lat_emission, lon_emission, nh3_emission_total_E

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

time_series, lat_emission, lon_emission, emission_mon_mean = monthly_mean_cal()
print(emission_mon_mean.shape,'!emission_mon_mean')

#compuation of annual emission
annual_emission = np.nansum(emission_mon_mean, axis=0) #in kg
print (annual_emission.shape, 'annual_emission.shape')

annual_emission[annual_emission <= 0.001] = np.nan


title_list = 'Annual Ammonia Emission + Vertical Profile boxes'
lat = lat_emission
lon = lon_emission

fig = plt.figure(facecolor='White',figsize=[11,14]);pad= 0.5; 
ax = plt.subplot(111);
plt.title(title_list, fontsize = 30, y=1.02)
colormap=cmap; colorbar_min=0;colorbar_max=0.20 ## change this accordingly
colormesh_1 = spatial_figure1(ax,annual_emission,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True)


ax.annotate('Wales',xy=(0.33,0.28), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
ax.annotate('S England',xy=(0.60,0.11), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)	
		
ax.annotate('E England',xy=(0.85,0.34), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)		

ax.annotate('NE England',xy=(0.65,0.54), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)			
		
ax.annotate('N Ireland',xy=(0.12,0.51), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)		

ax.annotate('Scotland',xy=(0.41,0.71), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points', color='black',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)		

cbar_ax = fig.add_axes([0.05, 0.05, 0.90, 0.02])  # position and size of the colorbar
char = fig.colorbar(colormesh_1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.tick_params(labelsize=25) 
cbar_ax.annotate(r'Annual NH$_3$ Emission (kg)',xy=(0.5,0.9), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=30)

plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'AnnualNH3_Emission_verticalProfileBoxes.png',bbox_inches='tight')   ##########	
plt.show()		
