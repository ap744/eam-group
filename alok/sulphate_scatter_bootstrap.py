# This code compare GEOS-Chem model and DEFRA sites sulphate 
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


#read UKEAP sulphate datasets here scratch_alok -> /scratch/uptrop/ap744
path='/scratch/uptrop/ap744/UKEAP_data/UKEAP_AcidGases_Aerosol/UKEAP_Particulate_Sulphate/'
sulphate_files=glob.glob(path + '28-UKA0*-2016_particulate_sulphate_*.csv')
print (sulphate_files)

# read csv file having DEFRA sites details
sites = pd.read_csv('/scratch/uptrop/ap744/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_AcidGases_Aerosol_sites_details.csv', encoding= 'unicode_escape')
#print (sites.head(10))
ID = sites["UK-AIR_ID"]
print (ID)

# site wise annual mean computation  
x = []
for f in sulphate_files:
	df = pd.read_csv(f,parse_dates=["Start Date", "End Date"])  
	print (df.head(5))
	print (len(sulphate_files))
	sitesA = sites.copy()
	#df['Measurement'].values[df['Measurement'] <=0.1] = np.nan

	#Annual Mean calculation
	mean_A= df["Measurement"].mean() # to compute annual mean
	print (mean_A, f[87:95])
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
	
	sitesA["sulphate_annual_mean"] = mean_A
	sitesA["sulphate_mam_mean"] = mean_mam
	sitesA["sulphate_jja_mean"] = mean_jja
	sitesA["sulphate_son_mean"] = mean_son
	sitesA["sulphate_djf_mean"] = mean_djf
	#print (sitesA.head(10))
	
	x.append(
	{
		'UK-AIR_ID':f[87:95],
		'sulphate_annual_mean':mean_A,
		'sulphate_mam_mean':mean_mam,
		'sulphate_jja_mean':mean_jja,
		'sulphate_son_mean':mean_son,
		'sulphate_djf_mean':mean_djf
		}
		)
	#print (x)
	
id_mean = pd.DataFrame(x)
#print (id_mean.head(3))

df_merge_col = pd.merge(sites, id_mean, on='UK-AIR_ID', how ='right')
print (df_merge_col.head(25))

#####export csv file having site wise annual mean information if needed 
#df_merge_col.to_csv(r'/home/a/ap744/scratch_alok/python_work/sulphate_annual_mean.csv')

#drop extra information from pandas dataframe
df_merge_colA = df_merge_col.drop(['S No','2016_Data'], axis=1)
print (df_merge_colA.head(5))

# change datatype to float to remove any further problems
df_merge_colA['Long'] = df_merge_colA['Long'].astype(float)
df_merge_colA['Lat'] = df_merge_colA['Lat'].astype(float)

#get sites information
sites_lon = df_merge_colA['Long']
sites_lat = df_merge_colA['Lat']
#getting annual mean data
sites_sulphate_AM = df_merge_colA['sulphate_annual_mean']

#seasonal mean data
sites_sulphate_mam = df_merge_colA['sulphate_mam_mean']
sites_sulphate_jja = df_merge_colA['sulphate_jja_mean']
sites_sulphate_son = df_merge_colA['sulphate_son_mean']
sites_sulphate_djf = df_merge_colA['sulphate_djf_mean']
sites_name = df_merge_colA['Site_Name']
print (sites_sulphate_AM, sites_name, sites_lat, sites_lon)

#####Reading GEOS-Chem files ################



#os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/AerosolMass/")
#Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))

#os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/SpeciesConc/")
#Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))


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
Aerosols = [xr.open_dataset(file) for file in Aerosols]

GC_surface_sulfate = [data['AerMassSO4'].isel(time=0,lev=0) for data in Aerosols]
#print (GC_surface_sulfate)

#Geos-Chem Annual Mean
GC_surface_sulfate_AM = sum(GC_surface_sulfate)/len(GC_surface_sulfate)
#print (GC_surface_sulfate_AM,'AnnualMean')
print (GC_surface_sulfate_AM.shape,'AnnualMean shape')

#Geos-Chem seasonal Mean
GC_surface_sulfate_mam = sum(GC_surface_sulfate[2:5])/len(GC_surface_sulfate[2:5])
#print (GC_surface_sulfate_mam.shape, 'MAM-shape')

GC_surface_sulfate_jja = sum(GC_surface_sulfate[5:8])/len(GC_surface_sulfate[5:8])
#print (GC_surface_sulfate_jja)

GC_surface_sulfate_son = sum(GC_surface_sulfate[8:11])/len(GC_surface_sulfate[8:11])
#print (GC_surface_sulfate_son)

GC_surface_sulfate_jf = sum(GC_surface_sulfate[0:2])/len(GC_surface_sulfate[0:2])
print (GC_surface_sulfate_jf, 'jf_shape')

GC_surface_sulfate_d = GC_surface_sulfate[11]
print (GC_surface_sulfate_d, 'd_shape')

#mean of JF and Dec using np.array --> creating problem in plotting
#GC_surface_sulfate_djf_a = np.array([GC_surface_sulfate_jf,GC_surface_sulfate_d])
#GC_surface_sulfate_djf = np.nanmean(GC_surface_sulfate_djf_a,axis=0)
#print (GC_surface_sulfate_djf, 'djf_shape')


GC_surface_sulfate_djf = (GC_surface_sulfate_d+GC_surface_sulfate_jf)/2
print (GC_surface_sulfate_djf, 'djf_shape')

#GEOS-Chem lat long information --Not working properly
#gc_lon = Aerosols[0]['lon']
#gc_lat = Aerosols[0]['lat']
#gc_lon,gc_lat = np.meshgrid(gc_lon,gc_lat)

# get GEOS-Chem lon and lat
gc_lon = GC_surface_sulfate_AM['lon']
gc_lat = GC_surface_sulfate_AM['lat']
print (len(gc_lon))
print (len(gc_lat))
print ((gc_lon))
print ((gc_lat))

# get number of sites from size of long and lat:
nsites=len(sites_lon)

# Define GEOS-Chem data obtained at same location as monitoring sites:
gc_data_sulphate_annual=np.zeros(nsites)

gc_data_sulphate_mam=np.zeros(nsites)
gc_data_sulphate_jja=np.zeros(nsites)
gc_data_sulphate_son=np.zeros(nsites)
gc_data_sulphate_djf=np.zeros(nsites)


#extract GEOS-Chem data using DEFRA sites lat long 
for w in range(len(sites_lat)):
	#print ((sites_lat[w],gc_lat))
	# lat and lon indices:
	lon_index = np.argmin(np.abs(np.subtract(sites_lon[w],gc_lon)))
	lat_index = np.argmin(np.abs(np.subtract(sites_lat[w],gc_lat)))

	#print (lon_index)
	#print (lat_index)
	gc_data_sulphate_annual[w] = GC_surface_sulfate_AM[lon_index, lat_index]
	gc_data_sulphate_mam[w] = GC_surface_sulfate_mam[lon_index, lat_index]
	gc_data_sulphate_jja[w] = GC_surface_sulfate_jja[lon_index, lat_index]
	gc_data_sulphate_son[w] = GC_surface_sulfate_son[lon_index, lat_index]
	gc_data_sulphate_djf[w] = GC_surface_sulfate_djf[lon_index, lat_index]

print (gc_data_sulphate_annual.shape)
print (sites_sulphate_AM.shape)

# quick scatter plot
#plt.plot(sites_sulphate_AM,gc_data_sulphate_annual,'o')
#plt.show()

# Compare DEFRA and GEOS-Chem:

#Normalized mean bias
nmb_Annual=100.*((np.nanmean(sites_sulphate_AM))- np.nanmean(gc_data_sulphate_annual))/np.nanmean(gc_data_sulphate_annual)
nmb_mam=100.*((np.nanmean(sites_sulphate_mam))- np.nanmean(gc_data_sulphate_mam))/np.nanmean(gc_data_sulphate_mam)
nmb_jja=100.*((np.nanmean(sites_sulphate_jja))- np.nanmean(gc_data_sulphate_jja))/np.nanmean(gc_data_sulphate_jja)
nmb_son=100.*((np.nanmean(sites_sulphate_son))- np.nanmean(gc_data_sulphate_son))/np.nanmean(gc_data_sulphate_son)
nmb_djf=100.*((np.nanmean(sites_sulphate_djf))- np.nanmean(gc_data_sulphate_djf))/np.nanmean(gc_data_sulphate_djf)
print(' DEFRA NMB_Annual= ', nmb_Annual)
print(' DEFRA NMB_mam = ', nmb_mam)
print(' DEFRA NMB_jja = ', nmb_jja)
print(' DEFRA NMB_son = ', nmb_son)
print(' DEFRA NMB_djf = ', nmb_djf)

nmb_Annual=100.*((np.nanmean(gc_data_sulphate_annual))- np.nanmean(sites_sulphate_AM))/np.nanmean(sites_sulphate_AM)
nmb_mam=100.*((np.nanmean(sites_sulphate_mam))- np.nanmean(gc_data_sulphate_mam))/np.nanmean(gc_data_sulphate_mam)
nmb_jja=100.*((np.nanmean(sites_sulphate_jja))- np.nanmean(gc_data_sulphate_jja))/np.nanmean(gc_data_sulphate_jja)
nmb_son=100.*((np.nanmean(sites_sulphate_son))- np.nanmean(gc_data_sulphate_son))/np.nanmean(gc_data_sulphate_son)
nmb_djf=100.*((np.nanmean(sites_sulphate_djf))- np.nanmean(gc_data_sulphate_djf))/np.nanmean(gc_data_sulphate_djf)
print(' DEFRA NMB_Annual= ', nmb_Annual)
print(' DEFRA NMB_mam = ', nmb_mam)
print(' DEFRA NMB_jja = ', nmb_jja)
print(' DEFRA NMB_son = ', nmb_son)
print(' DEFRA NMB_djf = ', nmb_djf)


#correlation
correlate_annual=stats.pearsonr(gc_data_sulphate_annual,sites_sulphate_AM)

# dropping nan values and compute correlation
nas_mam = np.logical_or(np.isnan(gc_data_sulphate_mam), np.isnan(sites_sulphate_mam))
correlate_mam = stats.pearsonr(gc_data_sulphate_mam[~nas_mam],sites_sulphate_mam[~nas_mam])

nas_jja = np.logical_or(np.isnan(gc_data_sulphate_jja), np.isnan(sites_sulphate_jja))
correlate_jja = stats.pearsonr(gc_data_sulphate_jja[~nas_jja],sites_sulphate_jja[~nas_jja])

nas_son = np.logical_or(np.isnan(gc_data_sulphate_son), np.isnan(sites_sulphate_son))
correlate_son = stats.pearsonr(gc_data_sulphate_son[~nas_son],sites_sulphate_son[~nas_son])

nas_djf = np.logical_or(np.isnan(gc_data_sulphate_djf), np.isnan(sites_sulphate_djf))
correlate_djf = stats.pearsonr(gc_data_sulphate_djf[~nas_djf],sites_sulphate_djf[~nas_djf])
print('Correlation = ',correlate_annual)

x11 = gc_data_sulphate_annual
y11 = sites_sulphate_AM
print (x11)
print (y11)
x11a=np.array([x11,y11])
x11b = x11a[~np.isnan(x11a).any(1)]
print (x11b.shape, 'X11b.shape')
print (x11b, 'X11b')
x1 = x11b[0]
y1 = x11b[1]
print (x1,'x1')
print (y1,'y1')


x22 = gc_data_sulphate_mam
y22 = sites_sulphate_mam
print (x22,'x22')
print (y22,'y22')
x22a=np.array([x22,y22])
print (x22a.shape, 'X22a.shape')
x22b = x22a[:,~np.isnan(x22a).any(axis=0)]
print (x22b, 'X22b')
x2 = x22b[0]
y2 = x22b[1]
print (x2, 'x2')
print (y2, 'y2')


x33 = gc_data_sulphate_jja
y33 = sites_sulphate_jja
x33a=np.array([x33,y33])
x33b = x33a[:,~np.isnan(x33a).any(axis=0)]
x3 = x33b[0]
y3 = x33b[1]

x44 = gc_data_sulphate_son
y44 = sites_sulphate_son
x44a=np.array([x44,y44])
x44b = x44a[:,~np.isnan(x44a).any(axis=0)]
x4 = x44b[0]
y4 = x44b[1]

x55 = gc_data_sulphate_djf
y55 = sites_sulphate_djf
x55a=np.array([x55,y55])
x55b = x55a[:,~np.isnan(x55a).any(axis=0)]
x5 = x55b[0]
y5 = x55b[1]



#Regression ~ using bootstrap
#ind=np.where((gc_data_sulphate_annual!=0)&(sites_sulphate_AM!=0))
#print (ind)
regres_annual=rma(x1,y1,(len(x1)),1000)
print('slope annual: ',regres_annual[0])   
print('Intercept annual: ',regres_annual[1])
print('slope error annual: ',regres_annual[2])   
print('Intercept error annual: ',regres_annual[3])


regres_mam=rma(x2,y2,(len(x2)),1000)
print('slope mam: ',regres_mam[0])   
print('Intercept mam: ',regres_mam[1])
print('slope error mam: ',regres_mam[2])   
print('Intercept error mam: ',regres_mam[3])

regres_jja=rma(x3,y3,(len(x3)),1000)
print('slope jja: ',regres_jja[0])   
print('Intercept jja: ',regres_jja[1])
print('slope error jja: ',regres_jja[2])   
print('Intercept error jja: ',regres_jja[3])


regres_son=rma(x4,y4,(len(x4)),1000)
print('slope son: ',regres_son[0])   
print('Intercept son: ',regres_son[1])
print('slope error son: ',regres_son[2])   
print('Intercept error son: ',regres_son[3])



regres_djf=rma(x5,y5,(len(x5)),1000)
print('slope djf: ',regres_djf[0])   
print('Intercept djf: ',regres_djf[1])
print('slope error djf: ',regres_djf[2])   
print('Intercept error djf: ',regres_djf[3])



#scatter plot 
title_list = 'DEFRA and GEOS-Chem Sulfate'
title_list1 = 'DEFRA & GEOS-Chem Sulphate Annual'
title_list2 = 'DEFRA & GEOS-Chem Sulphate MAM'
title_list3 = 'DEFRA & GEOS-Chem Sulphate JJA'
title_list4 = 'DEFRA & GEOS-Chem Sulphate SON'
title_list5 = 'DEFRA & GEOS-Chem Sulphate DJF'



#plotting scatter plot
fig1 = plt.figure(facecolor='White',figsize=[11,11]);pad= 1.1;

ax = plt.subplot(111);
plt.title(title_list1, fontsize = 20, y=1)

plt.plot(x1,y1, 'o', color='black',markersize=6)
xvals=np.arange(0,100,5)
yvals=regres_annual[1]+xvals*regres_annual[0]
plt.plot(xvals,yvals, '-', color = 'r')
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
plt.xlabel("GEOS-Chem sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
plt.legend(('Data','RegressionLine','1:1 Line',),fontsize = 16)
add2plt=("y = ({a:2.2f}±{c:2.2f})x + ({b:2.2f}±{d:2.2f})".\
		format(a=regres_annual[0],b=regres_annual[1],c=regres_annual[2],d=regres_annual[3]))
plt.text(1.1,1.7,add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R$^2$ = {a:6.2f}".format(a=correlate_annual[0]*correlate_annual[0]))
plt.text(1.1,1.6, add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'sulfate_GEOS-Chem_DEFRAscatter_annual_bootstrap.png',bbox_inches='tight')
#plt.show()

fig2 = plt.figure(facecolor='White',figsize=[15,18]);pad= 1.1;

ax = plt.subplot(221);
plt.title(title_list2, fontsize = 20, y=1)

plt.plot(x2,y2, 'o', color='black',markersize=6)
xvals=np.arange(0,100,5)
yvals=regres_mam[1]+xvals*regres_mam[0]
plt.plot(xvals,yvals, '-', color = 'r')
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
plt.xlabel("GEOS-Chem sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
plt.legend(('Data','RegressionLine','1:1 Line',),fontsize = 16)
add2plt=("y = ({a:2.2f}±{c:2.2f})x + ({b:2.2f}±{d:2.2f})".\
		format(a=regres_mam[0],b=regres_mam[1],c=regres_mam[2],d=regres_mam[3]))
plt.text(1.1,1.7,add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R$^2$ = {a:6.2f}".format(a=correlate_mam[0]*correlate_mam[0]))
plt.text(1.1,1.6, add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)
  
ax = plt.subplot(222);
plt.title(title_list3, fontsize = 20, y=1)

plt.plot(x3,y3, 'o', color='black',markersize=6)
xvals=np.arange(0,100,5)
yvals=regres_jja[1]+xvals*regres_jja[0]
plt.plot(xvals,yvals, '-', color = 'r')
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
plt.xlabel("GEOS-Chem sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
plt.legend(('Data','RegressionLine','1:1 Line',),fontsize = 16)
add2plt=("y = ({a:2.2f}±{c:2.2f})x + ({b:2.2f}±{d:2.2f})".\
		format(a=regres_jja[0],b=regres_jja[1],c=regres_jja[2],d=regres_jja[3]))
plt.text(1.1,1.7,add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R$^2$ = {a:6.2f}".format(a=correlate_jja[0]*correlate_jja[0]))
plt.text(1.1,1.6, add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)
  
ax = plt.subplot(223);
plt.title(title_list4, fontsize = 20, y=1)

plt.plot(x4,y4, 'o', color='black',markersize=6)
xvals=np.arange(0,100,5)
yvals=regres_son[1]+xvals*regres_son[0]
plt.plot(xvals,yvals, '-', color = 'r')
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
plt.xlabel("GEOS-Chem sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
plt.legend(('Data','RegressionLine','1:1 Line',),fontsize = 16)
add2plt=("y = ({a:2.2f}±{c:2.2f})x + ({b:2.2f}±{d:2.2f})".\
		format(a=regres_son[0],b=regres_son[1],c=regres_son[2],d=regres_son[3]))
plt.text(1.1,1.7,add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R$^2$ = {a:6.2f}".format(a=correlate_son[0]*correlate_son[0]))
plt.text(1.1,1.6, add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)


ax = plt.subplot(224);
plt.title(title_list5, fontsize = 20, y=1)

plt.plot(x5,y5, 'o', color='black',markersize=6)
xvals=np.arange(0,100,5)
yvals=regres_djf[1]+xvals*regres_djf[0]
plt.plot(xvals,yvals, '-', color = 'r')
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
plt.xlabel("GEOS-Chem sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_sulfate ($\mu$g m$^{-3}$)", fontsize = 20)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
plt.legend(('Data','RegressionLine','1:1 Line',),fontsize = 16)
add2plt=("y = ({a:2.2f}±{c:2.2f})x + ({b:2.2f}±{d:2.2f})".\
		format(a=regres_djf[0],b=regres_djf[1],c=regres_djf[2],d=regres_djf[3]))
plt.text(1.1,1.7,add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R$^2$ = {a:6.2f}".format(a=correlate_djf[0]*correlate_djf[0]))
plt.text(1.1,1.6, add2plt, fontsize=16,\
		ha='left', va='center')#, transform=ax.transAxes)
for axis in ['top','bottom','left','right']:
		ax.spines[axis].set_linewidth(2)


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'sulfate_GEOS-Chem_DEFRAscatter_seasonal_bootstrap.png',bbox_inches='tight')
plt.show()
