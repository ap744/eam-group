import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc4
import os
import glob
import iris
from iris.coords import DimCoord
from iris.cube import Cube
from mpl_toolkits.basemap import Basemap


IASI_PATH = '/scratch/uptrop/em440/for_Alok/iasi_ncdf/iasi_nh3_uk_oversampled_2008-2018_0.1_jul2020.nc'
NAEI_PATH = '/scratch/uptrop/em440/for_Alok/naei_nh3/NAEI_total_NH3_0.1x0.1_2016.nc'
GC_FOLDER_PATH = "/scratch/uptrop/em440/for_Alok/gc_ncdf/"
WORLD_PATH = '/scratch/uptrop/ap744/shapefiles/Shapfiles_india/World_shp/World'
FIGURE_DIR = '/home/j/jfr10/alok/figures'

SAT_UK_LAT_MIN_INDEX = 172
SAT_UK_LAT_MAX_INDEX = 279
SAT_UK_LON_MIN_INDEX = 50
SAT_UK_LON_MAX_INDEX = 176

IASI_UK_LAT_MIN_INDEX = 1
IASI_UK_LAT_MAX_INDEX = 108
IASI_UK_LON_MIN_INDEX = 1
IASI_UK_LON_MAX_INDEX = 127

NAEI_UK_LAT_MIN_INDEX = 7
NAEI_UK_LAT_MAX_INDEX = 114
NAEI_UK_LON_MIN_INDEX = 4
NAEI_UK_LON_MAX_INDEX = 130

REGRID_RES = 0.1

NA = 6.022e23  # molecules/mol
mNH3 = 17.0  # g(NO2)/mol
mair = 28.97  # g(air)/mol

ANALYSIS_YEAR = 2016

class BadBoundaryException(Exception):
	pass

def main():
	data_emission, gc_column = get_data_for_year(GC_FOLDER_PATH)

	sat_lats, sat_lons = get_lat_lon_scale(os.path.join(GC_FOLDER_PATH, "satellite_files", "ts_08_11.EU.20160101.nc"))
	data_emission, sat_lats, sat_lons = regrid(data_emission, sat_lats, sat_lons, REGRID_RES)
	gc_lats, gc_lons = get_lat_lon_scale(os.path.join(GC_FOLDER_PATH, "emissions", "HEMCO_diagnostics.201601010000.nc"))
	gc_column, gc_lats, gc_lons = regrid(gc_column, gc_lats, gc_lons, REGRID_RES)

	data_emission_uk = data_emission[SAT_UK_LAT_MIN_INDEX:SAT_UK_LAT_MAX_INDEX, SAT_UK_LON_MIN_INDEX:SAT_UK_LON_MAX_INDEX]
	gc_column_uk = gc_column[SAT_UK_LAT_MIN_INDEX:SAT_UK_LAT_MAX_INDEX, SAT_UK_LON_MIN_INDEX:SAT_UK_LON_MAX_INDEX]

	iasi_monthly_uk, iasi_lats, iasi_lons = read_variable_over_area(IASI_PATH, "iasi_nh3",
																	IASI_UK_LAT_MIN_INDEX, IASI_UK_LAT_MAX_INDEX,
																	IASI_UK_LON_MIN_INDEX, IASI_UK_LON_MAX_INDEX)
	
	naei_nh3, naei_lats, naei_lons = read_variable_over_area(NAEI_PATH, "NH3",
																	NAEI_UK_LAT_MIN_INDEX, NAEI_UK_LAT_MAX_INDEX,
																	NAEI_UK_LON_MIN_INDEX, NAEI_UK_LON_MAX_INDEX)
	naei_area, _, _ = read_variable_over_area(NAEI_PATH, "area",
																	NAEI_UK_LAT_MIN_INDEX, NAEI_UK_LAT_MAX_INDEX,
																	NAEI_UK_LON_MIN_INDEX, NAEI_UK_LON_MAX_INDEX)

	naei_monthly_uk = (naei_nh3 * naei_area )/1000   # kg/yr

	uk_mask = np.where(naei_monthly_uk >= 100, 1, 0)    # multiplicative mask

	annual_emission = np.nansum(data_emission_uk, axis=0)   # ALOK: SHOULD THIS BE JUST OVER THE UK?

	iasi_ratios, naei_ratios, differences = [], [], []
	for iasi_month, naei_month, gc_column_month, data_emission_month \
			in zip(iasi_monthly_uk, naei_monthly_uk, gc_column_uk, data_emission_uk):
		masked_iasi_ratio = calc_iasi_nh3(data_emission_month, gc_column_month, iasi_month) * uk_mask
		masked_naei_ratio = calc_naei_nh3(data_emission_month, naei_month, annual_emission) * uk_mask
		difference = masked_naei_ratio - masked_iasi_ratio
		iasi_ratios.append(masked_iasi_ratio)
		naei_ratios.append(masked_naei_ratio)
		differences.append(difference)

	iasi_plot = plt.figure()
	for month, iasi_ratio in enumerate(iasi_ratios):
		ax = plt.subplot(3, 4, month + 1)
		plot_dataset(iasi_ratio, ax, iasi_lats, iasi_lons)
	iasi_plot.savefig(os.path.join(FIGURE_DIR, 'IASI_derived_NH3emissionF11W.png'))
	
	naei_plot = plt.figure()
	for month, naei_ratio in enumerate(naei_ratios):
		ax = plt.subplot(3, 4, month + 1)
		plot_dataset(naei_ratio, ax, naei_lats, naei_lons)
	naei_plot.savefig(os.path.join(FIGURE_DIR, 'NAEI_NH3_emissionF11W.png'))

	diff_plot = plt.figure()
	for month, difference in enumerate(differences):
		ax = plt.subplot(3, 4, month + 1)
		plot_dataset(difference, ax, naei_lats, naei_lons)
	diff_plot.savefig(os.path.join(FIGURE_DIR, 'NAEI_IASI_difference_emissionF11W.png'))


def get_data_for_year(gc_folder):
	sat_stack_list, emission_stack_list = [], []
	for month in range(1,13):
		sat_mean, em_sum = get_data_for_month(gc_folder, month)
		sat_stack_list.append(sat_mean)
		emission_stack_list.append(em_sum)
	sat_year_means = np.dstack(sat_stack_list)
	em_year_sums = np.dstack(emission_stack_list)
	return sat_year_means, em_year_sums


def get_data_for_month(gc_folder, month):
	month_str = str.zfill(str(month), 2)
	satellite_glob = os.path.join(gc_folder, "satellite_files", "ts_08_11.EU.2016{}*.nc".format(month_str))
	satellite_file_list = sorted(glob.glob(satellite_glob))
	emissions_glob = os.path.join(gc_folder, "emissions", "HEMCO_diagnostics.2016{}*0000.nc".format(month_str))
	emissions_file_list = sorted(glob.glob(emissions_glob))
	sat_stack_list, emission_stack_list = [], []
	for sat_file_path, emission_file_path in zip(satellite_file_list, emissions_file_list):
		sat_stack_list.append(load_and_preproc_sat_data(sat_file_path))
		emission_stack_list.append(load_and_preproc_emission_data(emission_file_path))
	sat_stack = np.dstack(sat_stack_list)
	em_stack = np.dstack(emission_stack_list)
	sat_mean = np.nanmean(sat_stack, axis=2)
	em_sum = np.nansum(em_stack, axis=2)
	return sat_mean, em_sum


def load_and_preproc_sat_data(sat_file_path):
	ncf_sat = nc4.Dataset(sat_file_path, mode='r')
	nh3_gc_column = ncf_sat.variables['IJ-AVG-S__NH3'][:]
	airdensity_sat = ncf_sat.variables['TIME-SER__AIRDEN'][:]
	bxheight_sat = ncf_sat.variables['BXHGHT-S__BXHEIGHT'][:]
	ncf_sat.close()

	bxheight_sat = bxheight_sat*100
	airdensity_sat = airdensity_sat*bxheight_sat
	nh3_gc_column = (nh3_gc_column/1e9)*airdensity_sat
	return np.nansum(nh3_gc_column, axis=0)


def load_and_preproc_emission_data(em_file_path):
	ncf_em = nc4.Dataset(em_file_path, mode='r')
	nh3_emission_total = ncf_em.variables['EmisNH3_Total'][:]
	area_emissions = ncf_em.variables['AREA'][:]
	ncf_em.close()

	nh3_emission_total = nh3_emission_total[0, :, :, :]
	nh3_emission_total = np.nansum(nh3_emission_total, axis=0)/area_emissions
	return nh3_emission_total


def get_lat_lon_scale(data_path):
	data = nc4.Dataset(data_path, mode='r')
	try:
		lat = data.variables['LAT'][:]
		lon = data.variables['LON'][:]
	except KeyError:
		lat = data.variables['lat'][:]
		lon = data.variables['lon'][:]
	data.close()
	return lat, lon


def regrid(data, from_lat, from_lon, to_resolution):
	lat_min, lon_min = np.nanmin(from_lat), np.nanmin(from_lon)
	lat_max, lon_max = np.nanmax(from_lat), np.nanmax(from_lon)
	to_lat = np.arange(lat_min, lat_max, to_resolution)
	to_lon = np.arange(lon_min, lon_max, to_resolution)

	latitude = DimCoord(from_lat, standard_name='latitude', units='degrees')
	longitude = DimCoord(from_lon, standard_name='longitude', units='degrees')
	time = DimCoord(np.linspace(1, 12, 12), standard_name='time', units='month')
	cube1 = Cube(data, dim_coords_and_dims=[(time, 2), (latitude, 0), (longitude, 1)])
	regridded_data = cube1.interpolate([('latitude', to_lat), ('longitude', to_lon)], iris.analysis.Linear())
	reshaped_data = np.transpose(regridded_data.data[:], (2, 0, 1))
	return reshaped_data, regridded_data.coord('latitude').points[:], regridded_data.coord('longitude').points[:]


def read_variable_over_area(dataset_path, variable_id,
						lat_min_index, lat_max_index, lon_min_index, lon_max_index):
	lats, lons = get_lat_lon_scale(dataset_path)
	data_set = nc4.Dataset(dataset_path, mode='r')
	data = data_set.variables[variable_id][:]
	data_set.close()
	lat_view = lats[lat_min_index:lat_max_index]
	lon_view = lons[lon_min_index:lon_max_index]
	data_view = data[..., lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	if 0 in data_view.shape:
		print("Bad bounding box defined for {}".format(dataset_path))
		raise BadBoundaryException
	return data_view, lat_view, lon_view


def calc_iasi_nh3(data_emission, regridded_gc_column, iasi_column):
	gc_ratio = data_emission / regridded_gc_column
	iasi_nh3 = (gc_ratio * iasi_column)/1000   # I guess in kg?
	return iasi_nh3


def calc_naei_nh3(data_emission, naei_nh3_area, annual_emission):
	# NOTE TO SELF: WATCH OUT that naei_nh3_area is in kg/year
	scale_factor = data_emission/annual_emission
	naei_nh3 = (scale_factor * naei_nh3_area)/1000   # Again, assuming in kg?
	return naei_nh3


def discrete_cmap(N, base_cmap=None):
	"""Create an N-bin discrete colormap from the specified input map"""
	# Note that if base_cmap is a string or None, you can simply do
	#    return plt.cm.get_cmap(base_cmap, N)
	# The following works for string, None, or a colormap instance:
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0, 1, N))
	cmap_name = base.name + str(N)
	return base.from_list(cmap_name, color_list, N)


def spatial_figure(axs, data, lons, lats, colormap, colorbar_min, colorbar_max, tb_lef=True, tb_bot=True,
				   bad_data=False):  # c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
		lons and lats are 1-d array while data is 2-D array
		colorbar_min,colorbar_max specifies the minimum and maximum value you want to show with the hottest and coldest color respectively
		tb_lef and tb_bot specifies if you want to have axis labels, True fro yes
	output : a spatial map of the data
	"""

	lon_b = -9
	lon_e = 3  # lonW,lonE,

	lat_b = 49
	lat_e = 61  # latS,latN

	lon_bin = 4
	lat_bin = 4

	map = Basemap(lat_0=0, lon_0=0, llcrnrlon=lon_b, llcrnrlat=lat_b, urcrnrlon=lon_e, urcrnrlat=lat_e, ax=axs,
				  projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)

	if tb_lef:
		map.drawparallels(np.arange(round(lat_b, 0) - lat_bin, round(lat_e, 0) + lat_bin, lat_bin), labels=[1, 0, 0, 0],
						  linewidth=0.0, fontsize=24)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b, 0), round(lon_e, 0) + lon_bin, lon_bin), labels=[0, 0, 0, 1],
						  linewidth=0.0, fontsize=24)
	# Add Coastlines, States, and Country Boundaries
	# map.drawcoastlines(); map.drawcountries() #map.drawstates(); # draw border lines
	map.readshapefile(WORLD_PATH, 'World',
					  linewidth=2)  # to add a shapefile on a map
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20, colormap)  # use 20 color bins, this can be changed
	cmap.set_bad([1, 1, 1], alpha=1.0);
	if bad_data:
		# cmap.set_under('w');cmap.set_over('k')
		cmap.set_under('k');
		cmap.set_over('w')
	colormesh = map.pcolormesh(xi, yi, masked_obj, cmap=cmap, vmin=colorbar_min, vmax=colorbar_max, latlon=True)
	return colormesh


def plot_dataset(dataset, ax, title, lon_range, lat_range):
	pad = 1.1
	plt.title(title, fontsize=30, y=1)
	colormap = discrete_cmap();
	colorbar_min = 0;
	colorbar_max = 60  ## change this accordingly
	colormesh_1 = spatial_figure(ax, dataset, lon_range, lat_range, colormap, colorbar_min,
								 colorbar_max, tb_lef=True, tb_bot=True, bad_data=False)
	# ax.annotate('MAM',xy=(0.07,0.90), xytext=(0, pad),
	# xycoords='axes fraction', textcoords='offset points',
	# ha='center', va='bottom',rotation='horizontal',fontsize=30)
	ax.annotate('NAEI NH$_3$ = {0:.2f}'.format(title) + ' Gg', xy=(0.41, 0.001), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom', rotation='horizontal', fontsize=35, color='r')





if __name__ == "__main__":
	main()

