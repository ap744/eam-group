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
from bootstrap import rma
from scipy.stats import gaussian_kde

#Use todays date to plot files - good practice to save files name with date
Today_date=datetime.datetime.now().strftime("%Y%m%d")

###Different cmap options
# cmap = matplotlib.cm.get_cmap('brewer_RdBu_11')
# cmap = cm.jet
cmap = cm.rainbow
#cmap = cm.YlOrRd

def NH3():
	#read UKEAP ammonia datasets here scratch_alok -> /scratch/uptrop/ap744
	path='/scratch/uptrop/ap744/UKEAP_data/for_NH3andNH4ratio/NH3/'
	ammonia_files=glob.glob(path + '*.csv')
	#print (ammonia_files)
	print (len(ammonia_files))

	# read csv file having DEFRA sites details
	sites = pd.read_csv('/scratch/uptrop/ap744/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_NH3_sites_details.csv', encoding= 'unicode_escape')
	#print (sites.head(10))
	ID = sites["UK-AIR_ID"]
	#print (ID)

	# site wise annual mean computation  
	x = []
	for f in ammonia_files:
		df = pd.read_csv(f,parse_dates=["Start Date", "End Date"])  
		#print (df.head(5))
		#print (len(ammonia_files))
		sitesA = sites.copy()
		#df['Measurement'].values[df['Measurement'] <=0.1] = np.nan

		#Annual Mean calculation
		mean_A= df["Measurement"].mean() # to compute annual mean
		#print (mean_A, f[59:67])
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
		#print (mean_d, 'mean_d')
		
		
		jf_start = pd.to_datetime("01/01/2016")
		jf_end = pd.to_datetime("15/03/2016")
		jf_subset = df[(df["Start Date"] > jf_start) & (df["End Date"] < jf_end)]
		mean_jf = jf_subset["Measurement"].mean()
		#print (mean_jf, 'mean_jf')
		
		
		mean_djf_a  = np.array([mean_d, mean_jf])
		
		mean_djf = np.nanmean(mean_djf_a, axis=0)
		#print (mean_djf, 'mean_djf')
		
		sitesA["ammonia_annual_mean"] = mean_A
		sitesA["ammonia_mam_mean"] = mean_mam
		sitesA["ammonia_jja_mean"] = mean_jja
		sitesA["ammonia_son_mean"] = mean_son
		sitesA["ammonia_djf_mean"] = mean_djf
		#print (sitesA.head(50))
		
		x.append(
		{
			'UK-AIR_ID':f[59:67],
			'ammonia_annual_mean':mean_A,
			'ammonia_mam_mean':mean_mam,
			'ammonia_jja_mean':mean_jja,
			'ammonia_son_mean':mean_son,
			'ammonia_djf_mean':mean_djf
			}
			)
		#print (x)
		
	id_mean = pd.DataFrame(x)
	#print (id_mean.head(3))

	df_merge_col = pd.merge(sites, id_mean, on='UK-AIR_ID', how ='right')
	#print (df_merge_col.head(50))

	#####export csv file having site wise annual mean information if needed 
	#df_merge_col.to_csv(r'/home/a/ap744/scratch_alok/python_work/ammonia_annual_mean.csv')

	#drop extra information from pandas dataframe
	df_merge_colA = df_merge_col.drop(['S No'], axis=1)
	#print (df_merge_colA.head(50))

	# change datatype to float to remove any further problems
	df_merge_colA['Long'] = df_merge_colA['Long'].astype(float)
	df_merge_colA['Lat'] = df_merge_colA['Lat'].astype(float)

	#get sites information
	sites_lon = df_merge_colA['Long']
	sites_lat = df_merge_colA['Lat']
	#getting annual mean data
	sites_ammonia_AM = df_merge_colA['ammonia_annual_mean']

	#seasonal mean data
	sites_ammonia_mam = df_merge_colA['ammonia_mam_mean']
	sites_ammonia_jja = df_merge_colA['ammonia_jja_mean']
	sites_ammonia_son = df_merge_colA['ammonia_son_mean']
	sites_ammonia_djf = df_merge_colA['ammonia_djf_mean']
	sites_name = df_merge_colA['Site_Name']
	#print (sites_ammonia_AM, sites_name, sites_lat, sites_lon)


	#####Reading GEOS-Chem files ################
	os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/AerosolMass/2016/")
	Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))

	os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/SpeciesConc/2016/")
	Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))

	os.chdir("/scratch/uptrop/ap744/GEOS-Chem_outputs/")
	StateMet = sorted(glob.glob("GEOSChem.StateMet.2016*b.nc4"))

	Species = Species[:] 
	Aerosols = Aerosols[:]
	StateMet = StateMet[:]
	#print(Aerosols, Species, StateMet, sep = "\n")

	Species  = [xr.open_dataset(file) for file in Species]
	Aerosols = [xr.open_dataset(file) for file in Aerosols]
	StateMet = [xr.open_dataset(file) for file in StateMet]

	#ds = xr.open_mfdataset(StateMet)
	#monthly_data = ds.resample(time='m').mean()
	#print(monthly_data)
	#monthly_data_StateMet = StateMet.resample(freq = 'm', dim = 'time', how = 'mean')
	#print(monthly_data_StateMet)

	#ammonia sufrace layer
	GC_surface_ammonia = [data['SpeciesConc_NH3'].isel(time=0,lev=0) for data in Species]
	#print (GC_surface_ammonia)

	#Avogadro's number [mol-1]
	AVOGADRO = 6.022140857e+23

	# Typical molar mass of air [kg mol-1]
	MW_AIR = 28.9644e-3
	# convert unit for ammonia (dry mol/mol to ug/m3)
	surface_AIRDEN = [data['Met_AIRDEN'].isel(time=0,lev=0) for data in StateMet] #kg/m3
	surface_AIRNUMDEN_a = np.asarray(surface_AIRDEN)/MW_AIR #mol/m3
	surface_AIRNUMDEN_b = surface_AIRNUMDEN_a*AVOGADRO # unit molec air/m3
	surface_AIRNUMDEN = surface_AIRNUMDEN_b/1e6 #unit molec air/cm3

	surface_ammonia_mass  = [x*y*17/(6.022*1e11) for (x,y) in zip(GC_surface_ammonia,surface_AIRNUMDEN)]
	#print (surface_ammonia_mass)

	#Geos-Chem Annual Mean
	GC_surface_ammonia_AM = sum(surface_ammonia_mass)/len(surface_ammonia_mass)
	#print (GC_surface_ammonia_AM,'AnnualMean')
	#print (GC_surface_ammonia_AM.shape,'AnnualMean shape')

	#Geos-Chem seasonal Mean
	GC_surface_ammonia_mam = sum(surface_ammonia_mass[2:5])/len(surface_ammonia_mass[2:5])
	#print (GC_surface_ammonia_mam.shape, 'MAM-shape')

	GC_surface_ammonia_jja = sum(surface_ammonia_mass[5:8])/len(surface_ammonia_mass[5:8])
	#print (GC_surface_ammonia_jja)

	GC_surface_ammonia_son = sum(surface_ammonia_mass[8:11])/len(surface_ammonia_mass[8:11])
	#print (GC_surface_ammonia_son)

	GC_surface_ammonia_jf = sum(surface_ammonia_mass[0:2])/len(surface_ammonia_mass[0:2])
	#print (GC_surface_ammonia_jf, 'jf_shape')

	GC_surface_ammonia_d = surface_ammonia_mass[11]
	#print (GC_surface_ammonia_d, 'd_shape')

	#mean of JF and Dec using np.array --> creating problem in plotting
	#GC_surface_ammonia_djf_a = np.array([GC_surface_ammonia_jf,GC_surface_ammonia_d])
	#GC_surface_ammonia_djf = np.nanmean(GC_surface_ammonia_djf_a,axis=0)
	#print (GC_surface_ammonia_djf, 'djf_shape')


	GC_surface_ammonia_djf = (GC_surface_ammonia_d+GC_surface_ammonia_jf)/2
	#print (GC_surface_ammonia_djf, 'djf_shape')

	#GEOS-Chem lat long information --Not working properly
	#gc_lon = Aerosols[0]['lon']
	#gc_lat = Aerosols[0]['lat']
	#gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat)

	# get GEOS-Chem lon and lat
	gc_lon = GC_surface_ammonia_AM['lon']
	gc_lat = GC_surface_ammonia_AM['lat']
	#print (len(gc_lon))
	#print (len(gc_lat))
	#print ((gc_lon))
	#print ((gc_lat))

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

	print (gc_data_ammonia_annual.shape, 'gc_data_ammonia_annual')
	print (sites_ammonia_AM.shape, 'sites_ammonia_AM')

	return (GC_surface_ammonia_AM,GC_surface_ammonia_mam,GC_surface_ammonia_jja,GC_surface_ammonia_son,GC_surface_ammonia_djf,
		gc_data_ammonia_annual,gc_data_ammonia_mam,gc_data_ammonia_jja,gc_data_ammonia_son,gc_data_ammonia_djf,
		sites_ammonia_AM,sites_ammonia_mam,sites_ammonia_jja,sites_ammonia_son,sites_ammonia_djf)

#NH4
def NH4():
	#read UKEAP ammonium datasets here scratch_alok -> /scratch/uptrop/ap744
	path='/scratch/uptrop/ap744/UKEAP_data/for_NH3andNH4ratio/NH4/'
	ammonium_files=glob.glob(path + '*.csv')
	#print (ammonium_files)
	print (len(ammonium_files))
	# read csv file having DEFRA sites details
	sites = pd.read_csv('/scratch/uptrop/ap744/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_AcidGases_Aerosol_sites_details.csv', encoding= 'unicode_escape')
	#print (sites.head(10))
	ID = sites["UK-AIR_ID"]
	#print (ID)

	# site wise annual mean computation  
	x = []
	for f in ammonium_files:
		df = pd.read_csv(f,parse_dates=["Start Date", "End Date"])  
		#print (df.head(5))
		#print (len(ammonium_files))
		sitesA = sites.copy()
		#df['Measurement'].values[df['Measurement'] <=0.1] = np.nan

		#Annual Mean calculation
		mean_A= df["Measurement"].mean() # to compute annual mean
		#print (mean_A, f[59:67])
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
		#print (mean_d, 'mean_d')
		
		
		jf_start = pd.to_datetime("01/01/2016")
		jf_end = pd.to_datetime("15/03/2016")
		jf_subset = df[(df["Start Date"] > jf_start) & (df["End Date"] < jf_end)]
		mean_jf = jf_subset["Measurement"].mean()
		#print (mean_jf, 'mean_jf')
		
		
		mean_djf_a  = np.array([mean_d, mean_jf])
		
		mean_djf = np.nanmean(mean_djf_a, axis=0)
		#print (mean_djf, 'mean_djf')
		
		sitesA["ammonium_annual_mean"] = mean_A
		sitesA["ammonium_mam_mean"] = mean_mam
		sitesA["ammonium_jja_mean"] = mean_jja
		sitesA["ammonium_son_mean"] = mean_son
		sitesA["ammonium_djf_mean"] = mean_djf
		#print (sitesA.head(10))
		
		x.append(
		{
			'UK-AIR_ID':f[59:67],
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
	#print (df_merge_col.head(25))

	#####export csv file having site wise annual mean information if needed 
	#df_merge_col.to_csv(r'/home/a/ap744/scratch_alok/python_work/ammonium_annual_mean.csv')

	#drop extra information from pandas dataframe
	df_merge_colA = df_merge_col.drop(['S No','2016_Data'], axis=1)
	#print (df_merge_colA.head(5))

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
	#print (sites_ammonium_AM, sites_name, sites_lat, sites_lon)


	#####Reading GEOS-Chem files ################


	os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/AerosolMass/2016/")
	Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))

	os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/SpeciesConc/2016/")
	Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))

	os.chdir("/scratch/uptrop/ap744/GEOS-Chem_outputs/")
	StateMet = sorted(glob.glob("GEOSChem.StateMet.2016*b.nc4"))
	Species = Species[:] 
	Aerosols = Aerosols[:]
	#print(Aerosols, Species, sep = "\n")

	Species  = [xr.open_dataset(file) for file in Species]
	Aerosols = [xr.open_dataset(file) for file in Aerosols]

	GC_surface_ammonium = [data['AerMassSO4'].isel(time=0,lev=0) for data in Aerosols]
	#print (GC_surface_ammonium)

	#Geos-Chem Annual Mean
	GC_surface_ammonium_AM = sum(GC_surface_ammonium)/len(GC_surface_ammonium)
	#print (GC_surface_ammonium_AM,'AnnualMean')
	#print (GC_surface_ammonium_AM.shape,'AnnualMean shape')

	#Geos-Chem seasonal Mean
	GC_surface_ammonium_mam = sum(GC_surface_ammonium[2:5])/len(GC_surface_ammonium[2:5])
	#print (GC_surface_ammonium_mam.shape, 'MAM-shape')

	GC_surface_ammonium_jja = sum(GC_surface_ammonium[5:8])/len(GC_surface_ammonium[5:8])
	#print (GC_surface_ammonium_jja)

	GC_surface_ammonium_son = sum(GC_surface_ammonium[8:11])/len(GC_surface_ammonium[8:11])
	#print (GC_surface_ammonium_son)

	GC_surface_ammonium_jf = sum(GC_surface_ammonium[0:2])/len(GC_surface_ammonium[0:2])
	#print (GC_surface_ammonium_jf, 'jf_shape')

	GC_surface_ammonium_d = GC_surface_ammonium[11]
	#print (GC_surface_ammonium_d, 'd_shape')

	#mean of JF and Dec using np.array --> creating problem in plotting
	#GC_surface_ammonium_djf_a = np.array([GC_surface_ammonium_jf,GC_surface_ammonium_d])
	#GC_surface_ammonium_djf = np.nanmean(GC_surface_ammonium_djf_a,axis=0)
	#print (GC_surface_ammonium_djf, 'djf_shape')


	GC_surface_ammonium_djf = (GC_surface_ammonium_d+GC_surface_ammonium_jf)/2
	#print (GC_surface_ammonium_djf, 'djf_shape')

	#GEOS-Chem lat long information --Not working properly
	#gc_lon = Aerosols[0]['lon']
	#gc_lat = Aerosols[0]['lat']
	#gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat)

	# get GEOS-Chem lon and lat
	gc_lon = GC_surface_ammonium_AM['lon']
	gc_lat = GC_surface_ammonium_AM['lat']
	#print (len(gc_lon))
	#print (len(gc_lat))
	#print ((gc_lon))
	#print ((gc_lat))

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

	print (gc_data_ammonium_annual.shape, 'gc_data_ammonium_annual')
	print (sites_ammonium_AM.shape, 'sites_ammonium_AM')

	return(GC_surface_ammonium_AM,GC_surface_ammonium_mam,GC_surface_ammonium_jja,GC_surface_ammonium_son,GC_surface_ammonium_djf,
		gc_data_ammonium_annual,gc_data_ammonium_mam,gc_data_ammonium_jja,gc_data_ammonium_son,gc_data_ammonium_djf,
		sites_ammonium_AM,sites_ammonium_mam,sites_ammonium_jja,sites_ammonium_son,sites_ammonium_djf,sites_lon,sites_lat)

GC_surface_ammonia_AM,GC_surface_ammonia_mam,GC_surface_ammonia_jja,GC_surface_ammonia_son,GC_surface_ammonia_djf,gc_data_ammonia_annual,gc_data_ammonia_mam,gc_data_ammonia_jja,gc_data_ammonia_son,gc_data_ammonia_djf,sites_ammonia_AM,sites_ammonia_mam,sites_ammonia_jja, sites_ammonia_son,sites_ammonia_djf = NH3()
GC_surface_ammonium_AM,GC_surface_ammonium_mam,GC_surface_ammonium_jja,GC_surface_ammonium_son,GC_surface_ammonium_djf,gc_data_ammonium_annual,gc_data_ammonium_mam,gc_data_ammonium_jja,gc_data_ammonium_son,gc_data_ammonium_djf,sites_ammonium_AM,sites_ammonium_mam,sites_ammonium_jja,sites_ammonium_son,sites_ammonium_djf,sites_lon,sites_lat = NH4()

NH3_NH4_annual_sites = sites_ammonia_AM/sites_ammonium_AM
NH3_NH4_mam_sites = sites_ammonia_mam/sites_ammonium_mam
NH3_NH4_jja_sites = sites_ammonia_jja/sites_ammonium_jja
NH3_NH4_son_sites = sites_ammonia_son/sites_ammonium_son
NH3_NH4_djf_sites = sites_ammonia_djf/sites_ammonium_djf

NH3_NH4_annual_model = gc_data_ammonia_annual/gc_data_ammonium_annual
NH3_NH4_mam_model = gc_data_ammonia_mam/gc_data_ammonium_mam
NH3_NH4_jja_model = gc_data_ammonia_jja/gc_data_ammonium_jja
NH3_NH4_son_model = gc_data_ammonia_son/gc_data_ammonium_son
NH3_NH4_djf_model = gc_data_ammonia_djf/gc_data_ammonium_djf

#NHx calculation Geos-Chem Model for spatial plot
NHx_annual_GC = GC_surface_ammonia_AM+GC_surface_ammonium_AM
NHx_mam_GC = GC_surface_ammonia_mam+GC_surface_ammonium_mam
NHx_jja_GC = GC_surface_ammonia_jja+GC_surface_ammonium_jja
NHx_son_GC = GC_surface_ammonia_son+GC_surface_ammonium_son
NHx_djf_GC = GC_surface_ammonia_djf+GC_surface_ammonium_djf
#NHx calculation Geos-Chem Model at DEFRA sites
NHx_annual_model = gc_data_ammonia_annual + gc_data_ammonium_annual
NHx_mam_model = gc_data_ammonia_mam + gc_data_ammonium_mam
NHx_jja_model = gc_data_ammonia_jja + gc_data_ammonium_jja
NHx_son_model = gc_data_ammonia_son + gc_data_ammonium_son
NHx_djf_model = gc_data_ammonia_djf + gc_data_ammonium_djf
#NHx calculation for DEFRA sites
NHx_annual_sites = sites_ammonia_AM + sites_ammonium_AM
NHx_mam_sites = sites_ammonia_mam + sites_ammonium_mam
NHx_jja_sites = sites_ammonia_jja + sites_ammonium_jja
NHx_son_sites = sites_ammonia_son + sites_ammonium_son
NHx_djf_sites = sites_ammonia_djf + sites_ammonium_djf

#NH3/NHx calcualtion for Geso-Chem data
NH3NHxRatio_GC_annual = GC_surface_ammonia_AM/NHx_annual_GC
NH3NHxRatio_GC_mam = GC_surface_ammonia_mam/NHx_mam_GC
NH3NHxRatio_GC_jja = GC_surface_ammonia_jja/NHx_jja_GC
NH3NHxRatio_GC_son = GC_surface_ammonia_son/NHx_son_GC
NH3NHxRatio_GC_djf = GC_surface_ammonia_djf/NHx_djf_GC
#NH3/NHx calcualtion for Geso-Chem DEFRA sites data
NH3NHxRatio_model_annual = gc_data_ammonia_annual/NHx_annual_model
NH3NHxRatio_model_mam = gc_data_ammonia_mam/NHx_mam_model
NH3NHxRatio_model_jja = gc_data_ammonia_jja/NHx_jja_model
NH3NHxRatio_model_son = gc_data_ammonia_son/NHx_son_model
NH3NHxRatio_model_djf = gc_data_ammonia_djf/NHx_djf_model
#NH3/NHx calcualtion for DEFRA sites
NH3NHxRatio_sites_annual = sites_ammonia_AM/NHx_annual_sites
NH3NHxRatio_sites_mam = sites_ammonia_mam/NHx_mam_sites
NH3NHxRatio_sites_jja = sites_ammonia_jja/NHx_jja_sites
NH3NHxRatio_sites_son = sites_ammonia_son/NHx_son_sites
NH3NHxRatio_sites_djf = sites_ammonia_djf/NHx_djf_sites

nmb_Annual=100.*((np.nanmean(NH3NHxRatio_model_annual))- np.nanmean(NH3NHxRatio_sites_annual))/np.nanmean(NH3NHxRatio_sites_annual)
nmb_mam=100.*((np.nanmean(NH3NHxRatio_model_mam))- np.nanmean(NH3NHxRatio_sites_mam))/np.nanmean(NH3NHxRatio_sites_mam)
nmb_jja=100.*((np.nanmean(NH3NHxRatio_model_jja))- np.nanmean(NH3NHxRatio_sites_jja))/np.nanmean(NH3NHxRatio_sites_jja)
nmb_son=100.*((np.nanmean(NH3NHxRatio_model_son))- np.nanmean(NH3NHxRatio_sites_son))/np.nanmean(NH3NHxRatio_sites_son)
nmb_djf=100.*((np.nanmean(NH3NHxRatio_model_djf))- np.nanmean(NH3NHxRatio_sites_djf))/np.nanmean(NH3NHxRatio_sites_djf)
print(' DEFRA NMB_Annual= ', nmb_Annual)
print(' DEFRA NMB_mam = ', nmb_mam)
print(' DEFRA NMB_jja = ', nmb_jja)
print(' DEFRA NMB_son = ', nmb_son)
print(' DEFRA NMB_djf = ', nmb_djf)

#correlation
correlate_Annual=stats.pearsonr(NH3NHxRatio_model_annual,NH3NHxRatio_sites_annual)

# dropping nan values and compute correlation
nas_mam = np.logical_or(np.isnan(NH3NHxRatio_model_mam), np.isnan(NH3NHxRatio_sites_mam))
correlate_mam = stats.pearsonr(NH3NHxRatio_model_mam[~nas_mam],NH3NHxRatio_sites_mam[~nas_mam])

nas_jja = np.logical_or(np.isnan(NH3NHxRatio_model_jja), np.isnan(NH3NHxRatio_sites_jja))
correlate_jja = stats.pearsonr(NH3NHxRatio_model_jja[~nas_jja],NH3NHxRatio_sites_jja[~nas_jja])

nas_son = np.logical_or(np.isnan(NH3NHxRatio_model_son), np.isnan(NH3NHxRatio_sites_son))
correlate_son = stats.pearsonr(NH3NHxRatio_model_son[~nas_son],NH3NHxRatio_sites_son[~nas_son])

nas_djf = np.logical_or(np.isnan(NH3NHxRatio_model_djf), np.isnan(NH3NHxRatio_sites_djf))
correlate_djf = stats.pearsonr(NH3NHxRatio_model_djf[~nas_djf],NH3NHxRatio_sites_djf[~nas_djf])

# plotting spatial map model and DEFRA network 
os.chdir('/home/a/ap744/scratch_alok/shapefiles/GBP_shapefile')
Europe_shape = r'GBR_adm1.shp'
Europe_map = ShapelyFeature(Reader(Europe_shape).geometries(),
                               ccrs.PlateCarree(), edgecolor='black',facecolor='none')
print ('Shapefile_read')
title_list = 'DEFRA and GEOS-Chem ammonia (g)'
title_list1 = 'Spatial Map DEFRA and GEOS-Chem ammonia (g)'

#fig,ax = plt.subplots(2,1, figsize=(11,11))
fig1 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
# plt.suptitle(title_list, fontsize = 35, y=0.96)

ax = plt.subplot(231);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

NH3NHxRatio_GC_annual.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'NH3/NHx GEOS-Chem',
											'orientation':'horizontal'})

ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_annual,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_annual,
		cmap=cmap,s = 100,vmin = 0,vmax = 1)
		

ax.set_title('Geos-Chem & DEFRA NH3/NHx (Annual)')
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_Annual = {0:.2f}'.format(correlate_Annual[0]),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
ax.annotate('NMB Annual= {0:.2f}'.format(nmb_Annual),xy=(0.8,0.6), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
		
colorbar = plt.colorbar(PCM, ax=ax,label='NH3/NHx DEFRA',
                       orientation='vertical',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=10) 
colorbar.ax.xaxis.label.set_size(10)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'NH3_NHx_ratio_GeosChem_DEFRA_spatial_annual.png',bbox_inches='tight')
#plt.show()



fig2 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(232);
#plt.title(title_list1, fontsize = 30, y=1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

NH3NHxRatio_GC_mam.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'NH3/NHx GEOS-Chem',
											'orientation':'horizontal'})
											
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_mam,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_mam,
		cmap=cmap,s = 100,vmin = 0,vmax = 1)
		

ax.set_title('Geos-Chem & DEFRA NH3/NHx (mam)')
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_mam = {0:.2f}'.format(correlate_mam[0]),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
ax.annotate('NMB mam = {0:.2f}'.format(nmb_mam),xy=(0.8,0.6), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
		
colorbar = plt.colorbar(PCM, ax=ax,label='NH3/NHx DEFRA',
                       orientation='vertical',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=10) 
colorbar.ax.xaxis.label.set_size(10)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'NH3_NHx_ratio_GeosChem_DEFRA_spatial_mam.png',bbox_inches='tight')


#fig,ax = plt.subplots(2,1, figsize=(11,11))
fig3 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(233);
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

NH3NHxRatio_GC_jja.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'NH3/NHx GEOS-Chem',
											'orientation':'horizontal'})
											
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_jja,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_jja,
		cmap=cmap,s = 100,vmin = 0,vmax = 1)
		

ax.set_title('Geos-Chem & DEFRA NH3/NHx (jja)')
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_jja = {0:.2f}'.format(correlate_jja[0]),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
ax.annotate('NMB jja = {0:.2f}'.format(nmb_jja),xy=(0.8,0.6), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
		
colorbar = plt.colorbar(PCM, ax=ax,label='NH3/NHx DEFRA',
                       orientation='vertical',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=10) 
colorbar.ax.xaxis.label.set_size(10)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'NH3_NHx_ratio_GeosChem_DEFRA_spatial_jja.png',bbox_inches='tight')




fig4 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(234);
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

NH3NHxRatio_GC_son.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'NH3/NHx GEOS-Chem',
											'orientation':'horizontal'})
											
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_son,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_son,
		cmap=cmap,s = 100,vmin = 0,vmax = 1)
		

ax.set_title('Geos-Chem & DEFRA NH3/NHx (son)')
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_son = {0:.2f}'.format(correlate_son[0]),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
ax.annotate('NMB son = {0:.2f}'.format(nmb_son),xy=(0.8,0.6), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
		
colorbar = plt.colorbar(PCM, ax=ax,label='NH3/NHx DEFRA',
                       orientation='vertical',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=10) 
colorbar.ax.xaxis.label.set_size(10)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'NH3_NHx_ratio_GeosChem_DEFRA_spatial_son.png',bbox_inches='tight')


fig5 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;
ax = plt.subplot(235);
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(Europe_map)
ax.set_extent([-9, 3, 49, 61], crs=ccrs.PlateCarree()) # [lonW,lonE,latS,latN]

NH3NHxRatio_GC_djf.plot(ax=ax,cmap=cmap,vmin = 0,vmax =1,
								cbar_kwargs={'shrink': 0.5, 
											'pad' : 0.01,
											'label': 'NH3/NHx GEOS-Chem',
											'orientation':'horizontal'})
											
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_djf,
		facecolors='none',edgecolors='black',linewidths=5,s = 100)
ax.scatter(x=sites_lon, y=sites_lat,c=NH3NHxRatio_sites_djf,
		cmap=cmap,s = 100,vmin = 0,vmax = 1)
		

ax.set_title('Geos-Chem & DEFRA NH3/NHx (djf)')
PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes


ax.annotate('Correl_djf = {0:.2f}'.format(correlate_djf[0]),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
ax.annotate('NMB djf = {0:.2f}'.format(nmb_djf),xy=(0.8,0.6), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=15)
		
colorbar = plt.colorbar(PCM, ax=ax,label='NH3/NHx DEFRA',
                       orientation='vertical',shrink=0.5,pad=0.01)
colorbar.ax.tick_params(labelsize=10) 
colorbar.ax.xaxis.label.set_size(10)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'NH3_NHx_ratio_GeosChem_DEFRA_spatial_djf.png',bbox_inches='tight')


plt.show()

