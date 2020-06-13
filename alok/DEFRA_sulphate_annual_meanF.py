import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib inline
import os 
from sklearn.preprocessing import StandardScaler
import datetime
Today_date=datetime.datetime.now().strftime("%Y%m%d")
import matplotlib.cm as cm
import glob

path='/home/a/ap744/scratch_alok/UKEAP_data/UKEAP_AcidGases_Aerosol/UKEAP_Particulate_Sulphate/'
NO3files=glob.glob(path + '28-UKA0*-2016_particulate_sulphate_*.csv')

print (NO3files)

#sites = pd.read_csv('/home/a/ap744/scratch_alok/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_AcidGases_Aerosol_sites_details.csv', encoding= 'unicode_escape')
sites = pd.read_csv('/home/a/ap744/scratch_alok/UKEAP_data/DEFRA_UKEAP_sites_details/UKEAP_AcidGases_Aerosol_sites_details.csv', encoding= 'unicode_escape')
#print (sites.head(10))
ID = sites["UK-AIR_ID"]
print (ID)

x= []

for f in NO3files:
	df = pd.read_csv(f)  
	#print (df.head(5))
	print (len(NO3files))
	sitesA = sites.copy()
	#df['Measurement'].values[df['Measurement'] <=0.1] = np.nan

	mean_A= df["Measurement"].mean()
	print (mean_A, f[92:100])
	sitesA["sulphate_annual_mean"] = mean_A
	#print (sitesA.head(10))
	
	x.append(
	{
		'UK-AIR_ID':f[92:100],
		'sulphate_annual_mean':mean_A
		}
		)
	
	#print (x)
id_mean = pd.DataFrame(x)

#print (id_mean.head(3))

df_merge_col = pd.merge(sites, id_mean, on='UK-AIR_ID', how ='right')

print (df_merge_col.head(7))

df_merge_col.to_csv(r'/home/a/ap744/scratch_alok/python_work/sulphate_annual_mean.csv')
	
