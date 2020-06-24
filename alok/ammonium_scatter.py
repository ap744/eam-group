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


#####Reading GEOS-Chem files ################


#os.chdir("/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_eu_naei_iccw/AerosolMass/2016")
os.chdir("/scratch/uptrop/ap744/GEOS-Chem_outputs/")

Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))
Aerosols = sorted(glob.glob("GEOSChem.AerosolMass*nc4"))
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

GC_surface_ammonium_d = GC_surface_ammonium[0]
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
nmb_Annual=100.*((sites_ammonium_AM)- np.mean(gc_data_ammonium_annual))/(gc_data_ammonium_annual)
nmb_mam=100.*((sites_ammonium_mam)- np.mean(gc_data_ammonium_mam))/(gc_data_ammonium_mam)
nmb_jja=100.*((sites_ammonium_jja)- np.mean(gc_data_ammonium_jja))/(gc_data_ammonium_jja)
nmb_son=100.*((sites_ammonium_son)- np.mean(gc_data_ammonium_son))/(gc_data_ammonium_son)
nmb_djf=100.*((sites_ammonium_djf)- np.mean(gc_data_ammonium_djf))/(gc_data_ammonium_djf)
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

print('Correlation Annual = ',correlate_Annual)

#Regression ~ bootstrap is not working
#regres=rma(gc_data_ammonium,sites_ammonium_AM,1000)
#print('slope: ',regres[0])   
#print('Intercept: ',regres[1])
#print('slope error: ',regres[2])   
#print('Intercept error: ',regres[3])





#scatter plot 

title_list = 'DEFRA and GEOS-Chem Particulate ammonium'
title_list1 = 'DERFA & GEOS-Chem ammonium Annual'
title_list2 = 'DERFA & GEOS-Chem ammonium MAM'
title_list3 = 'DERFA & GEOS-Chem ammonium JJA'
title_list4 = 'DERFA & GEOS-Chem ammonium SON'
title_list5 = 'DERFA & GEOS-Chem ammonium DJF'

x11 = gc_data_ammonium_annual
y11 = sites_ammonium_AM
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


x22 = gc_data_ammonium_mam
y22 = sites_ammonium_mam
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


x33 = gc_data_ammonium_jja
y33 = sites_ammonium_jja
x33a=np.array([x33,y33])
x33b = x33a[:,~np.isnan(x33a).any(axis=0)]
x3 = x33b[0]
y3 = x33b[1]

x44 = gc_data_ammonium_son
y44 = sites_ammonium_son
x44a=np.array([x44,y44])
x44b = x44a[:,~np.isnan(x44a).any(axis=0)]
x4 = x44b[0]
y4 = x44b[1]

x55 = gc_data_ammonium_djf
y55 = sites_ammonium_djf
x55a=np.array([x55,y55])
x55b = x55a[:,~np.isnan(x55a).any(axis=0)]
x5 = x55b[0]
y5 = x55b[1]



# Stat Calculation
x1y1 = np.vstack([x1,y1])
z1 = gaussian_kde(x1y1)(x1y1)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1,y1)  
line1 = slope1*x1+intercept1

x2y2 = np.vstack([x2,y2])
z2 = gaussian_kde(x2y2)(x2y2)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2,y2)  
line2 = slope2*x2+intercept2

x3y3 = np.vstack([x3,y3])
z3 = gaussian_kde(x3y3)(x3y3)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x3,y3)  
line3 = slope3*x3+intercept3

x4y4 = np.vstack([x4,y4])
z4 = gaussian_kde(x4y4)(x4y4)
slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(x4,y4)  
line4 = slope4*x4+intercept4

x5y5 = np.vstack([x5,y5])
z5 = gaussian_kde(x5y5)(x5y5)
slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(x5,y5)  
line5 = slope5*x5+intercept5


#plotting scatter plot
fig1 = plt.figure(facecolor='White',figsize=[25,18]);pad= 1.1;



ax = plt.subplot(231);
plt.title(title_list1, fontsize = 20, y=1)
ax.scatter(x1, y1, c=z1, s=101, edgecolor='',cmap=cm.jet)
plt.plot(x1, line1, 'r-', linewidth=3)
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
# plt.axis('equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()

plt.xlabel("GEOS-Chem ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
# plt.legend(( 'RegressionLine R$^2$={0:.2f}'.format(r_valueA),'1:1 Line','Data'))
plt.legend(('RegressionLine' ,'1:1 Line','Data'),fontsize = 16)
# print slope1,round(slope1,2),str(round(intercept1,1))
ax.annotate('R$^2$ = {0:.2f}'.format(r_value1*r_value1),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
ax.annotate('y = {0:.2f}'.format(slope1)+'x + {0:.2f}'.format(intercept1),xy=(0.7,0.63), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
 
 

ax = plt.subplot(232);
plt.title(title_list2, fontsize = 20, y=1)
ax.scatter(x2, y2, c=z2, s=101, edgecolor='',cmap=cm.jet)
plt.plot(x2, line2, 'r-', linewidth=3)
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
# plt.axis('equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()

plt.xlabel("GEOS-Chem ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
# plt.legend(( 'RegressionLine R$^2$={0:.2f}'.format(r_valueA),'1:1 Line','Data'))
plt.legend(('RegressionLine' ,'1:1 Line','Data'),fontsize = 16)
# print slope1,round(slope1,2),str(round(intercept1,1))
ax.annotate('R$^2$ = {0:.2f}'.format(r_value2*r_value2),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
ax.annotate('y = {0:.2f}'.format(slope2)+'x + {0:.2f}'.format(intercept2),xy=(0.7,0.63), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
  
ax = plt.subplot(233);
plt.title(title_list3, fontsize = 20, y=1)
ax.scatter(x3, y3, c=z3, s=101, edgecolor='',cmap=cm.jet)
plt.plot(x3, line3, 'r-', linewidth=3)
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
# plt.axis('equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()

plt.xlabel("GEOS-Chem ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
# plt.legend(( 'RegressionLine R$^2$={0:.2f}'.format(r_valueA),'1:1 Line','Data'))
plt.legend(('RegressionLine' ,'1:1 Line','Data'),fontsize = 16)
# print slope1,round(slope1,2),str(round(intercept1,1))
ax.annotate('R$^2$ = {0:.2f}'.format(r_value3*r_value3),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
ax.annotate('y = {0:.2f}'.format(slope3)+'x + {0:.2f}'.format(intercept3),xy=(0.7,0.63), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)
  
ax = plt.subplot(234);
plt.title(title_list4, fontsize = 20, y=1)
ax.scatter(x4, y4, c=z4, s=101, edgecolor='',cmap=cm.jet)
plt.plot(x4, line4, 'r-', linewidth=3)
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
# plt.axis('equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()

plt.xlabel("GEOS-Chem ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
# plt.legend(( 'RegressionLine R$^2$={0:.2f}'.format(r_valueA),'1:1 Line','Data'))
plt.legend(('RegressionLine' ,'1:1 Line','Data'),fontsize = 16)
# print slope1,round(slope1,2),str(round(intercept1,1))
ax.annotate('R$^2$ = {0:.2f}'.format(r_value4*r_value4),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
ax.annotate('y = {0:.2f}'.format(slope4)+'x + {0:.2f}'.format(intercept4),xy=(0.7,0.63), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)


ax = plt.subplot(235);
plt.title(title_list5, fontsize = 20, y=1)
ax.scatter(x5, y5, c=z5, s=101, edgecolor='',cmap=cm.jet)
plt.plot(x5, line5, 'r-', linewidth=3)
lineStart = 0 
lineEnd = 2
plt.plot([lineStart, lineEnd], [lineStart, lineEnd], 'k-', color = 'c',linewidth=2)
plt.xlim(lineStart, lineEnd)
plt.ylim(lineStart, lineEnd)
# plt.axis('equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()

plt.xlabel("GEOS-Chem ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
plt.ylabel("DEFRA_ammonium ($\mu$g m$^{-3}$)", fontsize = 20)
ax.yaxis.get_offset_text().set_size(16)
plt.yticks(fontsize = 16)
ax.xaxis.get_offset_text().set_size(16)
plt.xticks(fontsize = 16)
# plt.legend(( 'RegressionLine R$^2$={0:.2f}'.format(r_valueA),'1:1 Line','Data'))
plt.legend(('RegressionLine' ,'1:1 Line','Data'),fontsize = 16)
# print slope1,round(slope1,2),str(round(intercept1,1))
ax.annotate('R$^2$ = {0:.2f}'.format(r_value5*r_value5),xy=(0.8,0.7), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
ax.annotate('y = {0:.2f}'.format(slope5)+'x + {0:.2f}'.format(intercept5),xy=(0.7,0.63), xytext=(0, pad),
		xycoords='axes fraction', textcoords='offset points',
		ha='center', va='bottom',rotation='horizontal',fontsize=20)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(3)


plt.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.97, wspace=0.20, hspace=0.05);
plt.savefig('/scratch/uptrop/ap744/python_work/'+Today_date+'ammonium_GEOS-Chem_DEFRAscatter.png',bbox_inches='tight')
plt.show()
