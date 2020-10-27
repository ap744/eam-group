# This code compare GEOS-Chem model and DEFRA sites ammonium 
# Please contact Alok Pandey ap744@leicester.ac.uk for any further clarifications or details

#import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from sklearn.preprocessing import StandardScaler
import datetime
import xarray as xr
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.cm as cm
import glob
from scipy import stats
#from bootstrap import rma
from scipy.stats import gaussian_kde

#Use todays date to plot files - good practice to save files name with date
Today_date=datetime.datetime.now().strftime("%Y%m%d")

###Different cmap options
# cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')
# cmap = cm.jet
cmap = cm.rainbow
#cmap = cm.YlOrRd

#read UKEAP ammonium datasets here scratch_alok -> /scratch/uptrop/ap744
path='/scratch/uptrop/ap744/UKEAP_data/UKEAP_NH3_particulate_ammonium/particulate_ammonium/'
ammonium_files=glob.glob(path + '27-UKA0*-2016_particulate_ammonium_*.csv')
print (ammonium_files)

# read csv file having DEFRA sites details
sites = pd.read_csv('/scratch/uptrop/ap744/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_AcidGases_Aerosol_sites_details.csv', encoding= 'unicode_escape')
#print (sites.head(10))
ID = sites["UK-AIR_ID"]
print (ID)

# site wise annual mean computation  
x = []
for f in ammonium_files:
	df = pd.read_csv(f,parse_dates=["Start Date", "End Date"])  
	print (df.head(5))
	print (len(ammonium_files))
	sitesA = sites.copy()
	#df['Measurement'].values[df['Measurement'] <=0.1] = np.nan

	#Annual Mean calculation
	mean_A= df["Measurement"].mean() # to compute annual mean
	print (mean_A, f[88:96])
		#MAM mean Calculation
	mam_start = pd.to_datetime("15/02/2016")
	mam_end = pd.to_datetime("15/06/2016")
	mam_subset = df[(df["Start Date"] > mam_start) & (df["End Date"] < mam_end)]
	mean_mam = mam_subset["Measurement"].mean()
	
	#JJA mean Calculation
	jja_start = pd.to_datetime("15/05/2016")
	jja_end = pd.to_datetime("15/09/2016")
	jja_subset = df[(df["Start Date"] > jja_start) & (df["End Date"] < jja_end)]
	mean_jja = jja_subset["Measurement"].mean()

	#SON mean Calculation
	son_start = pd.to_datetime("15/08/2016")
	son_end = pd.to_datetime("15/11/2016")
	son_subset = df[(df["Start Date"] > son_start) & (df["End Date"] < son_end)]
	mean_son = son_subset["Measurement"].mean()
	
	#DJF mean Calculation
	
	d_start = pd.to_datetime("15/11/2016")
	d_end = pd.to_datetime("31/12/2016")
	d_subset = df[(df["Start Date"] > d_start) & (df["End Date"] < d_end)]
	mean_d = d_subset["Measurement"].mean()
	print (mean_d, 'mean_d')
	
	
	jf_start = pd.to_datetime("01/01/2016")
	jf_end = pd.to_datetime("15/03/2016")
	jf_subset = df[(df["Start Date"] > jf_start) & (df["End Date"] < jf_end)]
	mean_jf = jf_subset["Measurement"].mean()
	print (mean_jf, 'mean_jf')
	
	
	mean_djf_a  = np.array([mean_d, mean_jf])
	
	mean_djf = np.nanmean(mean_djf_a, axis=0)
	print (mean_djf, 'mean_djf')
	
	sitesA["ammonium_annual_mean"] = mean_A
	sitesA["ammonium_mam_mean"] = mean_mam
	sitesA["ammonium_jja_mean"] = mean_jja
	sitesA["ammonium_son_mean"] = mean_son
	sitesA["ammonium_djf_mean"] = mean_djf
	#print (sitesA.head(10))
	
	x.append(
	{
		'UK-AIR_ID':f[88:96],
		'ammonium_annual_mean':mean_A,
		'ammonium_mam_mean':mean_mam,
		'ammonium_jja_mean':mean_jja,
		'ammonium_son_mean':mean_son,
		'ammonium_djf_mean':mean_djf
		}
		)
	#print (x)
	
id_mean = pd.DataFrame(x)
#print (id_mean.head(3))

df_merge_col = pd.merge(sites, id_mean, on='UK-AIR_ID', how ='right')
print (df_merge_col.head(25))

#####export csv file having site wise annual mean information if needed 
#df_merge_col.to_csv(r'/home/a/ap744/scratch_alok/python_work/ammonium_annual_mean.csv')

#drop extra information from pandas dataframe
df_merge_colA = df_merge_col.drop(['S No','2016_Data'], axis=1)
print (df_merge_colA.head(5))


###################################################################################
###########  Delete Data over Scotland           ##################################
###################################################################################
df_merge_colA.drop(df_merge_colA[df_merge_colA['Lat'] > 56].index, inplace = True) 
print(df_merge_colA.head(11)) 
df_merge_colA.reset_index(drop=True, inplace=True)
print(df_merge_colA.head(11)) 


# change datatype to float to remove any further problems
df_merge_colA['Long'] = df_merge_colA['Long'].astype(float)
df_merge_colA['Lat'] = df_merge_colA['Lat'].astype(float)

#get sites information
sites_lon = df_merge_colA['Long']
sites_lat = df_merge_colA['Lat']
#getting annual mean data
sites_ammonium_AM = df_merge_colA['ammonium_annual_mean']

#seasonal mean data
sites_ammonium_mam = df_merge_colA['ammonium_mam_mean']
sites_ammonium_jja = df_merge_colA['ammonium_jja_mean']
sites_ammonium_son = df_merge_colA['ammonium_son_mean']
sites_ammonium_djf = df_merge_colA['ammonium_djf_mean']
sites_name = df_merge_colA['Site_Name']
print (sites_ammonium_AM, sites_name, sites_lat, sites_lon)


"""#####Reading GEOS-Chem files ################
########################### 50% increase in NH3 Emission ##################################
os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_scale_nh3_emis/AerosolMass/2016/")
Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))

os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_scale_nh3_emis/SpeciesConc/2016/")
Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))

os.chdir("/scratch/uptrop/ap744/GEOS-Chem_outputs/")
StateMet = sorted(glob.glob("GEOSChem.StateMet.2016*b.nc4"))

Species = Species[:] 
Aerosols = Aerosols[:]
#print(Aerosols, Species, sep = "\n")

Species  = [xr.open_dataset(file) for file in Species]
Aerosols = [xr.open_dataset(file) for file in Aerosols]"""


##############  new to read files  #############
#####Reading GEOS-Chem files ################
path_AerosolMass_2 = "/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/AerosolMass/2016/"
########################### 50% increase in NH3 Emission ##################################
path_AerosolMass_50increase = "/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_scale_nh3_emis/AerosolMass/2016/"

os.chdir(path_AerosolMass_2)
Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))

Aerosols = Aerosols[:]
Aerosols = [xr.open_dataset(file) for file in Aerosols]



GC_surface_ammonium = [data['AerMassNH4'].isel(time=0,lev=0) for data in Aerosols]
#print (GC_surface_ammonium)

#Geos-Chem Annual Mean
GC_surface_ammonium_AM = sum(GC_surface_ammonium)/len(GC_surface_ammonium)
#print (GC_surface_ammonium_AM,'AnnualMean')
print (GC_surface_ammonium_AM.shape,'AnnualMean shape')

#Geos-Chem seasonal Mean
GC_surface_ammonium_mam = sum(GC_surface_ammonium[2:5])/len(GC_surface_ammonium[2:5])
#print (GC_surface_ammonium_mam.shape, 'MAM-shape')

GC_surface_ammonium_jja = sum(GC_surface_ammonium[5:8])/len(GC_surface_ammonium[5:8])
#print (GC_surface_ammonium_jja)

GC_surface_ammonium_son = sum(GC_surface_ammonium[8:11])/len(GC_surface_ammonium[8:11])
#print (GC_surface_ammonium_son)

GC_surface_ammonium_jf = sum(GC_surface_ammonium[0:2])/len(GC_surface_ammonium[0:2])
print (GC_surface_ammonium_jf, 'jf_shape')

GC_surface_ammonium_d = GC_surface_ammonium[11]
print (GC_surface_ammonium_d, 'd_shape')

#mean of JF and Dec using np.array --> creating problem in plotting
#GC_surface_ammonium_djf_a = np.array([GC_surface_ammonium_jf,GC_surface_ammonium_d])
#GC_surface_ammonium_djf = np.nanmean(GC_surface_ammonium_djf_a,axis=0)
#print (GC_surface_ammonium_djf, 'djf_shape')


GC_surface_ammonium_djf = (GC_surface_ammonium_d+GC_surface_ammonium_jf)/2
print (GC_surface_ammonium_djf, 'djf_shape')

#GEOS-Chem lat long information --Not working properly
#gc_lon = Aerosols[0]['lon']
#gc_lat = Aerosols[0]['lat']
#gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat)

# get GEOS-Chem lon and lat
gc_lon = GC_surface_ammonium_AM['lon']
gc_lat = GC_surface_ammonium_AM['lat']
print (len(gc_lon))
print (len(gc_lat))
print ((gc_lon))
print ((gc_lat))

# get number of sites from size of long and lat:
nsites=len(sites_lon)

# Define GEOS-Chem data obtained at same location as monitoring sites:
gc_data_ammonium_annual=np.zeros(nsites)

gc_data_ammonium_mam=np.zeros(nsites)
gc_data_ammonium_jja=np.zeros(nsites)
gc_data_ammonium_son=np.zeros(nsites)
gc_data_ammonium_djf=np.zeros(nsites)


#extract GEOS-Chem data using DEFRA sites lat long 
for w in range(len(sites_lat)):
	#print ((sites_lat[w],gc_lat))
	# lat and lon indices:
	lon_index = np.argmin(np.abs(np.subtract(sites_lon[w],gc_lon)))
	lat_index = np.argmin(np.abs(np.subtract(sites_lat[w],gc_lat)))

	#print (lon_index)
	#print (lat_index)
	gc_data_ammonium_annual[w] = GC_surface_ammonium_AM[lon_index, lat_index]
	gc_data_ammonium_mam[w] = GC_surface_ammonium_mam[lon_index, lat_index]
	gc_data_ammonium_jja[w] = GC_surface_ammonium_jja[lon_index, lat_index]
	gc_data_ammonium_son[w] = GC_surface_ammonium_son[lon_index, lat_index]
	gc_data_ammonium_djf[w] = GC_surface_ammonium_djf[lon_index, lat_index]

print (gc_data_ammonium_annual.shape)
print (sites_ammonium_AM.shape)

# quick scatter plot
#plt.plot(sites_ammonium_AM,gc_data_ammonium_annual,'o')
#plt.show()

# Compare DERFA and GEOS-Chem:

#Normalized mean bias
#nmb_Annual=100.*((np.nanmean(sites_ammonium_AM))- np.nanmean(gc_data_ammonium_annual))/np.nanmean(gc_data_ammonium_annual)
#nmb_mam=100.*((np.nanmean(sites_ammonium_mam))- np.nanmean(gc_data_ammonium_mam))/np.nanmean(gc_data_ammonium_mam)
#nmb_jja=100.*((np.nanmean(sites_ammonium_jja))- np.nanmean(gc_data_ammonium_jja))/np.nanmean(gc_data_ammonium_jja)
#nmb_son=100.*((np.nanmean(sites_ammonium_son))- np.nanmean(gc_data_ammonium_son))/np.nanmean(gc_data_ammonium_son)
#nmb_djf=100.*((np.nanmean(sites_ammonium_djf))- np.nanmean(gc_data_ammonium_djf))/np.nanmean(gc_data_ammonium_djf)

nmb_Annual=100.*((np.nanmean(gc_data_ammonium_annual))- np.nanmean(sites_ammonium_AM))/np.nanmean(sites_ammonium_AM)
nmb_mam=100.*((np.nanmean(gc_data_ammonium_mam))- np.nanmean(sites_ammonium_mam))/np.nanmean(sites_ammonium_mam)
nmb_jja=100.*((np.nanmean(gc_data_ammonium_jja))- np.nanmean(sites_ammonium_jja))/np.nanmean(sites_ammonium_jja)
nmb_son=100.*((np.nanmean(gc_data_ammonium_son))- np.nanmean(sites_ammonium_son))/np.nanmean(sites_ammonium_son)
nmb_djf=100.*((np.nanmean(gc_data_ammonium_djf))- np.nanmean(sites_ammonium_djf))/np.nanmean(sites_ammonium_djf)

print(' DEFRA NMB_Annual= ', nmb_Annual)
print(' DEFRA NMB_mam = ', nmb_mam)
print(' DEFRA NMB_jja = ', nmb_jja)
print(' DEFRA NMB_son = ', nmb_son)
print(' DEFRA NMB_djf = ', nmb_djf)
#correlation
correlate_Annual=stats.pearsonr(gc_data_ammonium_annual,sites_ammonium_AM)

# dropping nan values and compute correlation
nas_mam = np.logical_or(np.isnan(gc_data_ammonium_mam), np.isnan(sites_ammonium_mam))
correlate_mam = stats.pearsonr(gc_data_ammonium_mam[~nas_mam],sites_ammonium_mam[~nas_mam])

nas_jja = np.logical_or(np.isnan(gc_data_ammonium_jja), np.isnan(sites_ammonium_jja))
correlate_jja = stats.pearsonr(gc_data_ammonium_jja[~nas_jja],sites_ammonium_jja[~nas_jja])

nas_son = np.logical_or(np.isnan(gc_data_ammonium_son), np.isnan(sites_ammonium_son))
correlate_son = stats.pearsonr(gc_data_ammonium_son[~nas_son],sites_ammonium_son[~nas_son])

nas_djf = np.logical_or(np.isnan(gc_data_ammonium_djf), np.isnan(sites_ammonium_djf))
correlate_djf = stats.pearsonr(gc_data_ammonium_djf[~nas_djf],sites_ammonium_djf[~nas_djf])

print('Correlation = ',correlate_Annual)

#Regression ~ bootstrap 
#regres=rma(gc_data_ammonium,sites_ammonium_AM,1000)
#print('slope: ',regres[0])   
#print('Intercept: ',regres[1])
#print('slope error: ',regres[2])   
#print('Intercept error: ',regres[3])

# plotting spatial map model and DEFRA network 
os.chdir('/home/a/ap744/scratch_alok/shapefiles/GBP_shapefile')
Europe_shape = r'GBR_adm1.shp'
Europe_map = ShapelyFeature(Reader(Europe_shape).geometries(),
                               ccrs.PlateCarree(), edgecolor='black',facecolor='none')
print ('Shapefile_read')
title_list = 'DEFRA and GEOS-Chem Particulate ammonium'
title_list1 = 'Spatial Map DEFRA and GEOS-Chem Particulate ammonium'

#fig,ax = plt.subplots(2,1, figsize=(11,11))
fig1 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
# plt.suptitle(title_list, fontsize = 35, y=0.96)

ax = plt.subplot(231);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonium_AM.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1.5,
								cbar_kwargs={'shrink': 0.0, 
											'pad' : 0.09,
											'label': '',
											'orientation':'horizontal'})

ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_AM,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_AM,
		cmap=cmap,s = 100,vmin = 0,vmax = 1.5)
		
ax.set_title('DEFRA and GEOS-Chem ammonium (Annual)',fontsize=20)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_Annual = {0:.2f}'.format(correlate_Annual[0]),xy=(0.65,0.75), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
ax.annotate('NMB Annual= {0:.2f}'.format(nmb_Annual),xy=(0.65,0.85), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
		
colorbar = plt.colorbar(PCM, ax=ax,label='GEOS-Chem & DEFRA ammonium ($\mu$g m$^{-3}$)',
                        orientation='horizontal',shrink=0.5,pad=0.05)
colorbar.ax.tick_params(labelsize=20) 
colorbar.ax.xaxis.label.set_size(21)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAspatial_annualwithoutScotland_iccw.png',bbox_inches='tight')




fig2 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(232);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonium_mam.plot(ax=ax,cmap=cmap,vmin = 0,vmax =3,
								cbar_kwargs={'shrink': 0.0, 
											'pad' : 0.09,
											'label': '',
											'orientation':'horizontal'})
											
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_mam,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_mam,
		cmap=cmap,s = 100,vmin = 0,vmax = 3)
		
ax.set_title('DEFRA and GEOS-Chem ammonium (mam)',fontsize=20)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes



ax.annotate('Correl_mam = {0:.2f}'.format(correlate_mam[0]),xy=(0.65,0.75), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
ax.annotate('NMB mam= {0:.2f}'.format(nmb_mam),xy=(0.65,0.85), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
		
		
colorbar = plt.colorbar(PCM, ax=ax,label='DEFRA & GEOS-Chem ammonium ($\mu$g m$^{-3}$)',
                        orientation='horizontal',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=20) 
colorbar.ax.xaxis.label.set_size(21)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAspatial_mam_iccw.png',bbox_inches='tight')





#fig,ax = plt.subplots(2,1, figsize=(11,11))
fig3 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;

ax = plt.subplot(233);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonium_jja.plot(ax=ax,cmap=cmap,vmin = 0,vmax =3,
								cbar_kwargs={'shrink': 0.0, 
											'pad' : 0.09,
											'label': '',
											'orientation':'horizontal'})

ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_jja,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_jja,
		cmap=cmap,s = 100,vmin = 0,vmax = 3)
		
ax.set_title('DEFRA and GEOS-Chem ammonium (jja)',fontsize=20)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_jja = {0:.2f}'.format(correlate_jja[0]),xy=(0.65,0.75), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
ax.annotate('NMB jja= {0:.2f}'.format(nmb_jja),xy=(0.65,0.85), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
		
colorbar = plt.colorbar(PCM, ax=ax,label='DEFRA & GEOS-Chem ammonium ($\mu$g m$^{-3}$)',
                        orientation='horizontal',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=20) 
colorbar.ax.xaxis.label.set_size(21)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAspatial_jja_iccw.png',bbox_inches='tight')




fig4 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(234);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonium_son.plot(ax=ax,cmap=cmap,vmin = 0,vmax =3,
								cbar_kwargs={'shrink': 0.0, 
											'pad' : 0.09,
											'label': '',
											'orientation':'horizontal'})

ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_son,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_son,
		cmap=cmap,s = 100,vmin = 0,vmax = 3)
		
ax.set_title('DEFRA and GEOS-Chem ammonium (son)',fontsize=20)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_son = {0:.2f}'.format(correlate_son[0]),xy=(0.65,0.75), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
ax.annotate('NMB son = {0:.2f}'.format(nmb_son),xy=(0.65,0.85), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
		
colorbar = plt.colorbar(PCM, ax=ax,label='DEFRA & GEOS-Chem ammonium ($\mu$g m$^{-3}$)',
                        orientation='horizontal',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=20) 
colorbar.ax.xaxis.label.set_size(21)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAspatial_son_iccw.png',bbox_inches='tight')





fig5 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(235);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

GC_surface_ammonium_djf.plot(ax=ax,cmap=cmap,vmin = 0,vmax =3,
								cbar_kwargs={'shrink': 0.0, 
											'pad' : 0.09,
											'label': '',
											'orientation':'horizontal'})

ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_djf,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=sites_ammonium_djf,
		cmap=cmap,s = 100,vmin = 0,vmax = 3)
		
ax.set_title('DEFRA and GEOS-Chem ammonium (djf)',fontsize=20)
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_djf = {0:.2f}'.format(correlate_djf[0]),xy=(0.65,0.75), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
ax.annotate('NMB djf = {0:.2f}'.format(nmb_djf),xy=(0.65,0.85), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20,color ='w')
		
colorbar = plt.colorbar(PCM, ax=ax,label='DEFRA & GEOS-Chem ammonium ($\mu$g m$^{-3}$)',
                        orientation='horizontal',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=20) 
colorbar.ax.xaxis.label.set_size(21)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAspatial_djf_iccw.png',bbox_inches='tight')

plt.show()
