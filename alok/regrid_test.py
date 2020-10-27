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

def regridd_2d_L2H(lon_o,lat_o,data_o,lon_n,lat_n,lon_rpt,lat_rpt):
	def AreaWeight(lon1,lon2,lat1,lat2):
				'''calculate the earth radius in m2'''
				radius = 6371000;
				area = (np.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
				(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
				# print np.nansum(np.nansum(area,axis=1),axis=0)
				return area
	
	# Derive the original grids $area_o
	lon_i,lat_i=np.meshgrid(lon_o,lat_o)   # meshgrid both lon and lat into 2D fields
	lat_itv = lat_o[1] - lat_o[0];lon_itv = lon_o[1] - lon_o[0];  # This applies only for regular grids, otherwise change it to array
	area_rpt = AreaWeight(lon_i,lon_i+lon_itv,lat_i,lat_i+lat_itv)
	#GAo = AreaWeight(lon1_o,lon2_o,lat1_o,lat2_o) # Original grid   0.25*0.3125
	GAo = area_rpt
	
	Eo = data_o
	
	EoA =np.multiply(Eo ,GAo)	
	
	# Derive the target grids area $area_n
	lon_m,lat_m=np.meshgrid(lon_n,lat_n)   # meshgrid both lon and lat into 2D fields
	lat_res = lat_n[1] - lat_n[0];lon_res = lon_n[1] - lon_n[0];  # This applies only for regular grids, otherwise change it to array
	area_n = AreaWeight(lon_m,lon_m+lon_res,lat_m,lat_m+lat_res)
	#GAd = AreaWeight(lon1_d,lon2_d,lat1_d,lat2_d) # destination grid 0.1*0.1
	GAd = area_n
	
	Ef= EoA /(5*25)  # np.shape(Ef) == np.shape(Eo)# Ef has the same shape as the original data, you want to create an array with a shape which has the original rows*5, original columns*25
		
	data_rpt = np.repeat(np.repeat(Ef, repeats=5,axis=0),repeats=25,axis=1)
	
	No_box_lat=2
	No_box_lon=25
	
	mass = np.zeros((len(lat_n),len(lon_n)));
	for row in range(0,len(lat_n)):
					for column in range(0,len(lon_n)):            
									mass[row,column] =np.nansum(np.nansum(EoA[row*No_box_lat:(row+1)*No_box_lat,column*No_box_lon:(column+1)*No_box_lon],axis=1),axis=0)
	mass = np.divide(mass,GAd)
	return mass

emission_data_file = '/scratch/uptrop/ap744/python_work/20200724emission_monthly_GC.nc'
ncf_emission = nc4.Dataset(emission_data_file,mode='r')

lat_o = ncf_emission.variables['lat'][:]
lon_o = ncf_emission.variables['lon'][:]
nh3_emission = ncf_emission.variables['emission_nh3_GC'][:]  
print (lat_o.shape, 'lat_o.shape')
print (lon_o.shape, 'lon_o.shape')
print (nh3_emission.shape, 'nh3_emission_total.shape')

jan_nh3_emission = nh3_emission[0,:,:]

lat_n=np.arange(32.75,61.25,0.1)
lon_n= np.arange(-15,40,0.1)
print (lat_n.shape, 'lat_omi_regrid')
#print (lat_n); 
print (lon_n.shape, 'lon_omi_regrid')
#print (lon_n)
 
lat_rpt =100
lon_rpt =100

jan_regrid_emission = regridd_2d_L2H(lon_o,lat_o,jan_nh3_emission,lon_n,lat_n,lon_rpt,lat_rpt)		# Regrid
print (jan_regrid_emission.shape, 'jan_regrid_emission.shape')

print (np.nanmax(jan_nh3_emission), 'max_o')
print (np.nanmin(jan_nh3_emission), 'min_o')
print (np.nanmax(jan_regrid_emission), 'max_n')
print (np.nanmin(jan_regrid_emission), 'min_n')

title_list = 'NH3 Emission regridding comparison'					###########
title_list1 = 'Original'
title_list2 = 'Regridded'

fig = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1; 
plt.suptitle(title_list, fontsize = 25, y=1.001)

ax = plt.subplot(1,2,1);
plt.title(title_list1, fontsize = 30, y=1)
colormap=cmap; colorbar_min=1e-07;colorbar_max=1e-02## change this accordingly
colormesh_1 = spatial_figure(ax,jan_nh3_emission,lon_o,lat_o,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=True)
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
colormap=cmap; colorbar_min=0;colorbar_max=1## change this accordingly
colormesh_1 = spatial_figure(ax,jan_regrid_emission,lon_n,lat_n,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True,bad_data=True)
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
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'regrid.png',bbox_inches='tight')   ##########	
plt.show()
