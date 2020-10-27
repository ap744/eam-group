import os
import numpy as np
import pandas as pd
import xarray as xr

# build a "regridding" funtion


os.chdir('/scratch/uptrop/ap744/python_work/')
emission_file = '20200727emission_monthly_GC.nc'

emission = xr.open_dataset(emission_file)
print (emission)

emission1 = emission.to_dataframe()
print (emission1)

emission_Jan = emission.isel(time=[0])

print (emission_Jan)

# fisrt decide the target grids centres after regridding
out_lon = np.arange(-15,40+0.1,0.1)  # (lon_min,lon_max+resolution,lon_resolution)
out_lat = np.arange(32.75,61.25+0.1,0.1) # (lat_min,lat_max+resolution,lat_resolution)

def regrid_TROPOMI(data):
    data.insert(3, "lon_new", np.nan) # insert an empty column for the new lon 
    data.insert(4, "lat_new", np.nan) # insert an empty column for the new lat 
    for i in range(len(data)):        # for each point in your raw data, find its nearest point on the target grid
        data['lon_new'][i] = out_lon[np.argmin(abs(out_lon-data['lon'][i]))]
        data['lat_new'][i] = out_lat[np.argmin(abs(out_lat-data['lat'][i]))]
    data_regrid = data.groupby(['lon_new','lat_new'],as_index=False).mean() # get mean NH3 at the same grid centre
    data_regrid = data_regrid[['lon_new','lat_new','NH3']]
    return data_regrid
    
regrid_Jan = regrid_TROPOMI(emission1)
print (regrid_Jan)
