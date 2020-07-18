# This code compare GEOS-Chem model and DEFRA sites ammonia 
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


#####Reading GEOS-Chem files ################

#Avogadro's number [mol-1]
AVOGADRO = 6.022140857e+23

# Typical molar mass of air [kg mol-1]
MW_AIR = 28.9644e-3

os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/SpeciesConc/2016/")
Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))
#print (len(Species))
Species = Species[:] 
Species_1  = [xr.open_dataset(file) for file in Species]
#print (len(Species_1))
print (Species_1)
#ammonia sufrace layer
GC_surface_ammonia = [data['SpeciesConc_NH3'].isel(time=0,lev=0) for data in Species_1]
print (GC_surface_ammonia)




os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/StateMet/2016/")
StateMet = sorted(glob.glob("GEOSChem.StateMet*.nc4"))
print (len(StateMet))

StateMet = StateMet[:]
StateMet_1 = [xr.open_dataset(file) for file in StateMet]
#print ((StateMet_1[0]))
combined = xr.combine_nested(StateMet_1, concat_dim=("time"))
print (combined.indexes)
#monthly mean
StateMet_2=combined.groupby('time.month').mean()
#print (len(StateMet_2))
print (StateMet_2)
#StateMet_3 = list(StateMet_2.groupby('time'))
StateMet_3 = list(StateMet_2.groupby("month", squeeze=False))
print (StateMet_3)

# convert unit for ammonia (dry mol/mol to ug/m3)
surface_AIRDEN = [data['Met_AIRDEN'].isel(time=0,lev=0) for data in StateMet_3] #kg/m3

surface_AIRNUMDEN_a = np.asarray(surface_AIRDEN)/MW_AIR #mol/m3
surface_AIRNUMDEN_b = surface_AIRNUMDEN_a*AVOGADRO # unit molec air/m3
surface_AIRNUMDEN = surface_AIRNUMDEN_b/1e6 #unit molec air/cm3

surface_ammonia_mass  = [x*y*17/(6.022*1e11) for (x,y) in zip(GC_surface_ammonia,surface_AIRNUMDEN)]
print (surface_ammonia_mass)

#Geos-Chem Annual Mean
GC_surface_ammonia_AM = sum(surface_ammonia_mass)/len(surface_ammonia_mass)
#print (GC_surface_ammonia_AM,'AnnualMean')
print (GC_surface_ammonia_AM.shape,'AnnualMean shape')

#Geos-Chem seasonal Mean
GC_surface_ammonia_mam = sum(surface_ammonia_mass[2:5])/len(surface_ammonia_mass[2:5])
#print (GC_surface_ammonia_mam.shape, 'MAM-shape')

GC_surface_ammonia_jja = sum(surface_ammonia_mass[5:8])/len(surface_ammonia_mass[5:8])
#print (GC_surface_ammonia_jja)

GC_surface_ammonia_son = sum(surface_ammonia_mass[8:11])/len(surface_ammonia_mass[8:11])
#print (GC_surface_ammonia_son)

GC_surface_ammonia_jf = sum(surface_ammonia_mass[0:2])/len(surface_ammonia_mass[0:2])
print (GC_surface_ammonia_jf, 'jf_shape')

GC_surface_ammonia_d = surface_ammonia_mass[11]
print (GC_surface_ammonia_d, 'd_shape')

#mean of JF and Dec using np.array --> creating problem in plotting
#GC_surface_ammonia_djf_a = np.array([GC_surface_ammonia_jf,GC_surface_ammonia_d])
#GC_surface_ammonia_djf = np.nanmean(GC_surface_ammonia_djf_a,axis=0)
#print (GC_surface_ammonia_djf, 'djf_shape')


GC_surface_ammonia_djf = (GC_surface_ammonia_d+GC_surface_ammonia_jf)/2
print (GC_surface_ammonia_djf, 'djf_shape')

#GEOS-Chem lat long information --Not working properly
#gc_lon = Aerosols[0]['lon']
#gc_lat = Aerosols[0]['lat']
#gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat)

# get GEOS-Chem lon and lat
gc_lon = GC_surface_ammonia_AM['lon']
gc_lat = GC_surface_ammonia_AM['lat']
print (len(gc_lon))
print (len(gc_lat))
print ((gc_lon))
print ((gc_lat))

# get number of sites from size of long and lat:
nsites=len(sites_lon)

# Define GEOS-Chem data obtained at same location as monitoring sites:
gc_data_ammonia_annual=np.zeros(nsites)

gc_data_ammonia_mam=np.zeros(nsites)
gc_data_ammonia_jja=np.zeros(nsites)
gc_data_ammonia_son=np.zeros(nsites)
gc_data_ammonia_djf=np.zeros(nsites)


#extract GEOS-Chem data using DEFRA sites lat long 
for w in range(len(sites_lat)):
	#print ((sites_lat[w],gc_lat))
	# lat and lon indices:
	lon_index = np.argmin(np.abs(np.subtract(sites_lon[w],gc_lon)))
	lat_index = np.argmin(np.abs(np.subtract(sites_lat[w],gc_lat)))

	#print (lon_index)
	#print (lat_index)
	gc_data_ammonia_annual[w] = GC_surface_ammonia_AM[lon_index, lat_index]
	gc_data_ammonia_mam[w] = GC_surface_ammonia_mam[lon_index, lat_index]
	gc_data_ammonia_jja[w] = GC_surface_ammonia_jja[lon_index, lat_index]
	gc_data_ammonia_son[w] = GC_surface_ammonia_son[lon_index, lat_index]
	gc_data_ammonia_djf[w] = GC_surface_ammonia_djf[lon_index, lat_index]

print (gc_data_ammonia_annual.shape)
print (sites_ammonia_AM.shape)
