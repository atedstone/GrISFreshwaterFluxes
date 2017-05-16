import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt

import georaster

pstere = {'init':'epsg:3413'}

## Load ice discharge datasets

## King/Howat/OSU monthly dataset

# Create column headings 
king_cols = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/monthly_discharge.csv',
	skiprows=[1, 3], nrows=2)

kcols = ['Year', 'Month']
skip = 0
for k, i in king_cols.iteritems():
	if skip < 2:
		skip += 1
		continue

	if k.find('Unnamed') == -1:
		name = k.replace(' ', '').lower().strip()
		name_col = name + '_discharge'
	else:
		name_col = name + '_std'

	kcols.append(name_col)

# Use columns to load dataset
king = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/monthly_discharge.csv',
	skiprows=4, names=kcols, parse_dates={'date':['Year', 'Month']}, index_col='date')


# King glaciers locations
king_locs = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/monthly_discharge.csv',
	skiprows=0, nrows=1)
king_locs.columns = kcols
lats = king_locs.filter(like='discharge')
lons = king_locs.filter(like='std')
lons.columns = lats.columns
king_locs = pd.concat([lats.T, lons.T], axis=1)
king_locs.columns = ['latitude', 'longitude']
geometry = [Point(xy) for xy in zip(king_locs.longitude, king_locs.latitude)]
king_geo = gpd.GeoDataFrame({'king_name':king_locs.index}, crs={'init':'epsg:4326'}, geometry=geometry)
#king_geo.geometry = king_geo.buffer(0.1)
king_geo.to_file('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/king_locs.shp')
king_pstere = king_geo.to_crs(pstere)

## Enderlin annual dataset
# Cut the various average rows off the bottom of the dataset
enderlin_raw = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/GrIS_D-timeseries.txt', encoding='mac_roman', 
	delim_whitespace=True, index_col='Glacier_Name', nrows=178)

enderlin = enderlin_raw.drop(enderlin_raw.columns[0], axis=1)
enderlin = enderlin.drop(enderlin_raw.columns[1], axis=1)
enderlin = enderlin.T
enderlin.index = pd.date_range('2000-01-01', '2012-01-01', freq='1AS')
enderlin_cols = [c + '_discharge' for c in enderlin.columns]
enderlin.columns = enderlin_cols


enderlin_locs = pd.concat([enderlin_raw[enderlin_raw.columns[0]], enderlin_raw[enderlin_raw.columns[1]]], axis=1)
enderlin_locs.columns = ['longitude', 'latitude']
enderlin_locs.index = enderlin_cols
geometry = [Point(xy) for xy in zip(enderlin_locs.longitude, enderlin_locs.latitude)]
enderlin_geo = gpd.GeoDataFrame({'enderlin_name':enderlin_locs.index}, crs={'init':'epsg:4326'}, geometry=geometry)
enderlin_geo.to_file('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/enderlin_locs.shp')
enderlin_pstere = enderlin_geo.to_crs(pstere)

## Find equivalent labels between Enderlin and King datasets (metres - using pstere coords)
king_equiv = []
king_equiv_ix = []
for ix, row in king_pstere.iterrows():
	dists = enderlin_pstere.distance(row.geometry)
	if dists.min() < 15000:
		ix_min = dists.argmin()
		equiv_label = enderlin_pstere.loc[ix_min].enderlin_name
		print('King: %s, Enderlin: %s' % (row.king_name, equiv_label))
		king_equiv.append(row.king_name)
		king_equiv_ix.append(equiv_label)
	else:
		king_equiv.append(np.nan)
		king_equiv_ix.append(np.nan)

# At this point I had to do some manual re-mapping of a few labels...
# ...hence the manually-specified list of names below.
enderlin_names_for_king = ['79fjorden_discharge',
	'hayes_discharge',
	'zachariaeisstrom_discharge',
	'alison_discharge',
	'daugaard_jensen_discharge',
	'helheim_discharge',
	'humboldt_discharge',
	'ikertivaq_south_discharge',
	'jakobshavn_discharge',
	'kangerdlugssuaq_discharge',
	'kangerdlugssup_discharge',
	'kangilerngata_discharge',
	'koge_bugt_discharge',
	'kongoscar_discharge',
	'petermann_discharge',
	'rink_discharge',
	'store_discharge',
	'tingmjarmiut_discharge',
	'torsukatat_discharge',
	'ukassorssuaq_discharge',
	'upernavik_central_discharge',
	'upernavik_north_discharge',
	'upernavik_northwest_discharge',
	'upernavik_south_discharge'
	]
# It seems impossible to match nordre* glaciers, so drop them here (and dropped in list above already)
king = king.drop('nordrenorth_discharge', axis=1)
king = king.drop('nordresouth_discharge', axis=1)



# Convert King's monthly rate to monthly flux then calculate total ice-sheet-wide annual flux
king = king.assign(year_length=np.where(king.index.is_leap_year, 366, 365))
king = king.assign(month_length=king.index.daysinmonth)
monthly_rate = king.filter(like='discharge')
monthly_flux = monthly_rate.apply(lambda x: (x / king.year_length.values) * king.month_length.values)
king_annual_total_flux = monthly_flux.resample('1AS').sum().T.sum()



## Compare Enderlin and King
king_annual_glacier_flux = monthly_flux.resample('1AS').sum()
for nk, ne in zip(king.filter(like='discharge').columns, enderlin_names_for_king):
	print(nk, ne)	
	plt.figure()
	plt.title(ne)
	plt.plot(enderlin.index, enderlin[ne], 'r')
	plt.plot(king_annual_glacier_flux.index, king_annual_glacier_flux[nk], 'b')
	plt.show()

# Discharge record from two glaciers does not match - ukassorssuaq and torsukatat
# Nevertheless, continue for the moment
# Enderlin total for King glaciers:
enderlin_annual_subset = enderlin.filter(items=enderlin_names_for_king).T.sum()
print(enderlin_annual_subset)
print(king_annual_total_flux)


## Calculate % contribution of each month to annual flux...
# percentages make glaciers comparable to one another
# first convert to monthly perc of annual each year at each glacier
monthly_perc = monthly_flux.groupby(monthly_flux.index.year).apply(lambda x: (100 / x.sum()) * x)
monthly_perc.groupby(monthly_perc.index.month).mean()
## the monthly percentage is actually very similar every month on average, i.e. c. 8%/year.




### Rignot data set

rignot_raw = pd.read_excel('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/Bamber+co-ords_reformatted.xlsx', index_col='rignot_name')
ix = [pd.datetime(y, 1, 1) for y in rignot_raw.columns[2:]]
rignot = rignot_raw.drop(rignot_raw.columns[0], axis=1)
rignot = rignot.drop(rignot_raw.columns[1], axis=1)
rignot = rignot.T
rignot.index = ix

rignot_glaciers = []
for c in rignot.columns:
	c = c.lower().replace(' ', '') + '_discharge'
	rignot_glaciers.append(c)

rignot.columns = rignot_glaciers

# Coordinates in spreadsheet are i,j IDL format, bamber grid, 2.5 km spacing
# Use lookup matrices to identify pstere coordinates
bamber_grid_y, bamber_grid_x = np.mgrid[-3400000:-600000:2500, -800000:700000:2500]

# For debugging use - can write out a raster grid to check it
# grid_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
# trans = (bamber_grid_x[0,0], (bamber_grid_x[0,1]-bamber_grid_x[0,0]), 0, bamber_grid_y[0,0], 0, (bamber_grid_y[1,0]-bamber_grid_y[0,0]))
# georaster.simple_write_geotiff(
# 	'bamber_grid_y.tif',
# 	bamber_grid_y,
# 	trans,
# 	proj4=grid_proj.srs,
# 	dtype=gdal.GDT_Float32
# 	)

rignot_x = []
rignot_y = []
for x, y in zip(rignot_raw['x'], rignot_raw['y']):
	rignot_y.append(bamber_grid_y[y, x])
	rignot_x.append(bamber_grid_x[y, x])

rignot_locs = pd.DataFrame({'x':rignot_x, 'y':rignot_y})
# stored as kms, so convert to metres

rignot_locs.index = rignot_glaciers
geometry = [Point(xy) for xy in zip(rignot_locs.x, rignot_locs.y)]
rignot_proj4 = '+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
rignot_bamber = gpd.GeoDataFrame({'rignot_name':rignot_locs.index}, crs=rignot_proj4, geometry=geometry)
rignot_geo = rignot_bamber.to_crs({'init':'epsg:4326'})
rignot_geo.to_file('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/rignot_locs.shp')
rignot_pstere = rignot_geo.to_crs(pstere)

rignot_equiv = []
rignot_equiv_ix = []
for ix, row in rignot_pstere.iterrows():
	dists = enderlin_pstere.distance(row.geometry)
	if dists.min() < 50000:
		ix_min = dists.argmin()
		equiv_label = enderlin_pstere.loc[ix_min].enderlin_name
		print('Rignot: %s, Enderlin: %s' % (row.rignot_name, equiv_label))
		rignot_equiv.append(row.rignot_name)
		rignot_equiv_ix.append(equiv_label)
	else:
		rignot_equiv.append(np.nan)
		rignot_equiv_ix.append(np.nan)


geodf_all = enderlin_pstere.merge(pd.DataFrame({'enderlin_name':rignot_equiv}, index=rignot_equiv), on='enderlin_name')

enderlin_master = pd.DataFrame({'enderlin_dummy':1}, index=enderlin_pstere.enderlin_name)
geodf_all = pd.merge(
	enderlin_master, 
	pd.DataFrame({'rignot_name':rignot_equiv}, index=rignot_equiv_ix), 
	how='left', 
	left_index=True,
	right_index=True)

geodf_all = pd.merge(
	geodf_all, 
	pd.DataFrame({'king_name':king_equiv}, index=king_equiv_ix), 
	how='left', 
	left_index=True,
	right_index=True)

geodf_all = geodf_all.drop('enderlin_dummy', axis=1)
geodf_all = geodf_all.drop('geometry', axis=1)

geodf_all.to_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/glacier_names_lookup.csv')





## Resolve individual outlet glaciers to the basins used in Rignot's 2008 analysis
rignot_basins = gpd.read_file('/home/at15963/Dropbox/work/papers/bamber_fwf/JLB_Analysis_2015/combined.shp')
# remove the 'basin' prefix
basins_pstere = rignot_basins.to_crs(pstere)
basins_pstere.columns = ['GRIDCODE', 'area', 'geometry', 'basin']
enderlin_names_basins = gpd.sjoin(enderlin_pstere, basins_pstere, how='left', op='within')

# For the outlets that are not within a basin, resolve them to the nearest one
enderlin_names_basins = enderlin_names_basins.assign(dist=0)
unalloc = enderlin_names_basins[enderlin_names_basins.basin.isnull()]
for ix, row in unalloc.iterrows():
	dists = basins_pstere.distance(row.geometry)
	ix_min = dists.argmin()
	nearest_basin = basins_pstere.loc[ix_min, 'basin']
	enderlin_names_basins.loc[ix, 'basin'] = nearest_basin	
	enderlin_names_basins.loc[ix, 'dist'] = dists.min()


## Calculate basin-by-basin fluxes
grouped = enderlin_names_basins.groupby('basin')
store = {}
for basin, glaciers in grouped:
	q = enderlin.filter(items=glaciers.enderlin_name).sum(axis=1)
	store[basin] = q
basin_q_enderlin = pd.DataFrame(store)
basin_q_enderlin = basin_q_enderlin.assign(basin38=basin_q_enderlin.filter(like='basin38_').sum(axis=1))
basin_q_enderlin = basin_q_enderlin.drop(labels=basin_q_enderlin.filter(like='basin38_').columns.tolist(), axis=1)

## Rignot by basin (rather than by glacier name, but it's the same thing really)
rignot_raw2 = pd.read_excel('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/Bamber+co-ords_updatedbasins.xlsx', index_col='rignot_name')
# ix = [pd.datetime(y, 1, 1) for y in rignot_raw.columns[2:]]
# rignot = rignot_raw.drop(rignot_raw.columns[0], axis=1)
# rignot = rignot.drop(rignot_raw.columns[1], axis=1)
# rignot = rignot.T
# rignot.index = ix

rignot_names_basins = pd.Series(rignot_raw2['basin'])
rignot_names_basins.index = rignot_glaciers
rignot_names_basins = rignot_names_basins[rignot_names_basins.notnull()]
# Convert from float to string with basin prefix
rignot_names_basins = pd.Series(['basin{:.0f}'.format(n) for n in rignot_names_basins], index=rignot_names_basins.index)
# Un-intuitive: if ends with b then returns true, otherwise null, so want to retain null values
#rignot_names_basins = rignot_names_basins[rignot_names_basins.str.endswith('b').isnull()]
grouped = rignot_names_basins.groupby(rignot_names_basins)
store = {}
for basin, glaciers in grouped:
	q = rignot.filter(items=glaciers.index).sum(axis=1)
	store[basin] = q
basin_q_rignot = pd.DataFrame(store)

# facet plot of basin-by-basin comparisons
figure()
n = 1
for b in basin_q_rignot.columns:
	subplot(8, 6, n)
	title(b)
	try:
		plot(basin_q_rignot.index, basin_q_rignot[b], 'blue')
		plot(basin_q_enderlin.index, basin_q_enderlin[b], 'red')
		xlim('2000-01-01', '2015-12-31')
	except KeyError:
		pass
	n += 1


# Use Enderlin data to work out % contribution of each glacier to basin's outflow
grouped = enderlin_names_basins.groupby('basin')
store = {}
for basin, glaciers in grouped:
	# Calculate mean annual basin-wide Q
	q_basin = enderlin.filter(items=glaciers.enderlin_name).sum(axis=1).mean()
	# Calculate mean annual Q per glacier
	q_glaciers = enderlin.filter(items=glaciers.enderlin_name).mean(axis=0)
	# Calculate %contribution of each glacier to basin Q
	qp = (100 / q_basin) * q_glaciers
	# Store it
	store[basin] = qp
glacier_basin_contrib = pd.DataFrame(store)

# Split Rignot per-basin data out to Enderlin-defined outlets
# Produces a 'per-glacier' time series like Enderlin's
for ix, row in rignot.iterrows():


"""
Next steps here:
finish up association of Rignot glaciers with Enderlin
 	- see QGIS file
compare Rignot and Enderlin, remove magnitude difference
for each glacier, create complete time series of Rignot-Enderlin-King
for coordinates of each glacier use Enderlin (at least as far as possible)
correlate annual spatial mean of entire time series with runoff
Use correlation to extend discharge time series
	- for temporal (monthly) distn use the monthly percentage calculated from the King dataset
	- not yet sure how to to distribute total annual pre-92 discharge spatially?
"""	





# joined = gpd.sjoin(enderlin_geo, king_geo, how='right', op='within')

# len(joined)



# # Compare the King and Enderlin datasets for common glaciers
# king_annual_glacier_flux = monthly_flux.resample('1AS').sum()
# for ne, nk in zip(joined.enderlin_name, joined.king_name):
# 	plt.figure()
# 	plt.title(ne)
# 	plt.plot(enderlin.index, enderlin[ne], 'r')
# 	plt.plot(king_annual_glacier_flux.index, king_annual_glacier_flux[nk], 'b')	