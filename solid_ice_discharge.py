import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt

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
king_geo.geometry = king_geo.buffer(0.1)




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



## Find equivalent labels between Enderlin and King datasets
for ix, row in king_geo.iterrows():
	dists = enderlin_geo.distance(row.geometry)
	ix_min = dists.argmin()
	equiv_label = enderlin_geo.loc[ix_min].enderlin_name
	print('King: %s, Enderlin: %s' % (row.king_name, equiv_label))

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
glacier_annual_flux = monthly_flux.resample('1AS').sum()
monthly_perc = 
monthly_flux.groupby(monthly_flux.index.month)




# joined = gpd.sjoin(enderlin_geo, king_geo, how='right', op='within')

# len(joined)



# # Compare the King and Enderlin datasets for common glaciers
# king_annual_glacier_flux = monthly_flux.resample('1AS').sum()
# for ne, nk in zip(joined.enderlin_name, joined.king_name):
# 	plt.figure()
# 	plt.title(ne)
# 	plt.plot(enderlin.index, enderlin[ne], 'r')
# 	plt.plot(king_annual_glacier_flux.index, king_annual_glacier_flux[nk], 'b')	