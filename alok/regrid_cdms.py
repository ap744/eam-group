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
import cdms2,regrid2
import cdat_info

"""emission_data_file = '/scratch/uptrop/ap744/python_work/20200724emission_monthly_GC.nc'
ncf_emission = nc4.Dataset(emission_data_file,mode='r')
lat_o = ncf_emission.variables['lat'][:]
lon_o = ncf_emission.variables['lon'][:]
nh3_emission = ncf_emission.variables['emission_nh3_GC'][:]  
print (lat_o, lat_o.shape, 'lat_o.shape')
print (lon_o, lon_o.shape, 'lon_o.shape')
print (nh3_emission.shape, 'nh3_emission_total.shape')

lat_uk = lat_o[66:]
print (lat_uk, lat_uk.shape, 'lat_uk.shape')

lon_uk = lon_o[16:64]
print (lon_uk, lon_uk.shape, 'lon_uk.shape')

UK_nh3 = nh3_emission[:,66:,16:64]
print (UK_nh3.shape)"""



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




f1 = cdms2.open('/scratch/uptrop/ap744/python_work/20200724emission_monthly_GC.nc')
f2 = cdms2.open('/scratch/uptrop/em440/for_Alok/naei_nh3/NAEI_total_NH3_0.1x0.1_2016.nc')

nh3_gc = f1('emission_nh3_GC')
print (nh3_gc.shape, 'nh3_GC.shape')

nh3_gc_a = nh3_gc[:,66:,16:64]
print (nh3_gc_a.shape, 'nh3_gc_a.shape')

lat_o = nh3_gc_a.getLatitude()[:]
lon_o = nh3_gc_a.getLongitude()[:]
print (lat_o, lon_o, 'lat-lon_o shape')

jan_nh3_emission = nh3_gc_a[4]
print (jan_nh3_emission.shape, 'originaldata shape')

#plt.imshow(nh3_gc[4])
#plt.show()

nh3_NAEI = f2('NH3')
print (nh3_NAEI.shape, 'nh3_NAEI.shape')
#print("Latitude values:",nh3_NAEI.getLatitude()[:-6])
#print("Longidute values:",nh3_NAEI.getLongitude()[:-2])

nh3_NAEI_a = nh3_gc[:,0:-6,0:-2]
print (nh3_NAEI_a.shape, 'nh3_NAEI_shape after slicing.shape')

outgrid = nh3_NAEI_a.getGrid()
print (outgrid, 'outgrid.shape')

nh3_gc_new =  nh3_gc_a.regrid(outgrid)
print (nh3_gc_new.shape, 'regridded.shape')
#plt.imshow(nh3_gc_new[4])
#plt.show()

lat_n = nh3_gc_new.getLatitude()[:]
lon_n = nh3_gc_new.getLongitude()[:]
print (lat_n, lon_n, 'lat lon new shape')
jan_regrid_emission = nh3_gc_new[4]
print (jan_regrid_emission.shape, 'regridded_emission_shape')


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
colormap=cmap; colorbar_min=1e-07;colorbar_max=1e-02## change this accordingly
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
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'regrid11.png',bbox_inches='tight')   ##########	
plt.show()




