"""
Generates a monthly solid ice discharge time series for the Greenland Ice 
Sheet, 1958 to 2016.

Uses two datasets: Rignot, King. 
- Rignot basins are scaled by King fluxes and used to provide data for 1958, 1964, 1992-1999 period.
- King is 2000-2016 'gold standard' dataset by-glacier, monthly

Andrew Tedstone (a.j.tedstone@bristol.ac.uk), May 2017 (version 1) then December 2017 (revised)
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import statsmodels.api as sm
from scipy import spatial
import georaster

from matplotlib import rcParams
rcParams['font.size'] = 8
rcParams['font.sans-serif'] = 'Arial'
rcParams['figure.titlesize'] = 8

pstere = {'init':'epsg:3413'}
plot_figs = True

## Processing path
# Laptop Ubuntu VM
folder_path = '/media/sf_C_DRIVE/Users/at15963/Dropbox/work/papers/bamber_fwf/'
# Desktop
# folder_path = '~/Dropbox/work/papers/bamber_fwf/'
### Load ice discharge datasets

## King/Howat/OSU monthly dataset
# First create/tidy up the column headings
king_cols = pd.read_csv(folder_path + 'monthlyGrISdischarge.csv',
	skiprows=[1, 3], nrows=2)

kcols = ['Year', 'Month']
skip = 0
for k, i in king_cols.iteritems():
	if skip < 2:
		skip += 1
		continue

	if k.find('[]') == -1:
		name = k.replace(' ', '').lower().strip().strip("'")
		name_col = name + '_discharge'
	else:
		name_col = name + '_std'

	kcols.append(name_col)

# Now use columns to load dataset
king = pd.read_csv(folder_path + 'monthlyGrISdischarge.csv',
	skiprows=4, names=kcols, parse_dates={'date':['Year', 'Month']}, index_col='date')


# King glaciers locations
king_locs = pd.read_csv(folder_path + 'monthlyGrISdischarge.csv',
	skiprows=2, nrows=1, names=kcols)
#king_locs.columns = kcols
lats = king_locs.filter(like='discharge')
lons = king_locs.filter(like='std')
lons.columns = lats.columns
king_locs = pd.concat([lats.T, lons.T], axis=1)
king_locs.columns = ['latitude', 'longitude']
geometry = [Point(xy) for xy in zip(king_locs.longitude, king_locs.latitude)]
king_geo = gpd.GeoDataFrame({'king_name':king_locs.index}, crs={'init':'epsg:4326'}, geometry=geometry)
king_geo.to_file(folder_path + 'outputs_Nov2017/king_locs.shp')
king_pstere = king_geo.to_crs(pstere)




# Retain only discharge, get rid of std
king = king.filter(like='_discharge')


## Convert King's monthly rate to monthly flux then calculate total ice-sheet-wide annual flux
king = king.assign(year_length=np.where(king.index.is_leap_year, 366, 365))
king = king.assign(month_length=king.index.daysinmonth)
monthly_rate = king.filter(like='discharge')
monthly_flux = monthly_rate.apply(lambda x: (x / king.year_length.values) * king.month_length.values)
# Annual flux ice-sheet-wide
king_annual_total_flux = monthly_flux.resample('1AS').sum().T.sum()
# Annual flux per flacier
king_annual_glacier_flux = monthly_flux.resample('1AS').sum()


print(king_annual_total_flux)


## Calculate % contribution of each month to annual flux...
# percentages make glaciers comparable to one another
# first convert to monthly perc of annual each year at each glacier - i.e. still a time series
monthly_perc = monthly_flux.groupby(monthly_flux.index.year).apply(lambda x: (100 / x.sum()) * x)
#monthly_perc.groupby(monthly_perc.index.month).mean()
## the monthly percentage is actually very similar every month on average, i.e. c. 8%/year.



## Load Rignot's drainage basins
rignot_basins = gpd.read_file(folder_path + 'JLB_Analysis_2015/combined.shp')
# remove the 'basin' prefix
basins_pstere = rignot_basins.to_crs(pstere)
basins_pstere.columns = ['GRIDCODE', 'area', 'geometry', 'basin']


## Resolve King individual outlet glaciers into their Rignot basins
# First do a spatial join, based on glacier point *within* basin poly
king_names_basins = gpd.sjoin(king_pstere, basins_pstere, how='left', op='within')

## Next, for the outlets that are not within a basin, resolve them to the nearest one
# Add a distance attribute
king_names_basins = king_names_basins.assign(dist=0)
unalloc = king_names_basins[king_names_basins.basin.isnull()]
for ix, row in unalloc.iterrows():
	dists = basins_pstere.distance(row.geometry)
	ix_min = dists.argmin()
	nearest_basin = basins_pstere.loc[ix_min, 'basin']
	king_names_basins.loc[ix, 'basin'] = nearest_basin	
	king_names_basins.loc[ix, 'dist'] = dists.min()


## For mass continuity with Rignot, move some sub-basins into bigger basins...
# allocate basins below Helheim into Helheim basin
king_names_basins.loc[king_names_basins.basin.str.contains('basin46'), 'basin'] = 'basin11'
king_names_basins.loc[king_names_basins.basin.str.contains('basin47'), 'basin'] = 'basin11'
# ikertivaq north, pamiataq, unnamed into ikertivaq (as sub-basin is on flow divide)
king_names_basins.loc[king_names_basins.basin.str.contains('basin48'), 'basin'] = 'basin12'
king_names_basins.loc[king_names_basins.basin.str.contains('basin49'), 'basin'] = 'basin12'
# Allocate all glaciers in basin38_* basins to basin38
king_names_basins.loc[king_names_basins.basin.str.contains('basin38_'), 'basin'] = 'basin38'


## Aggregate King measurements to basin outputs
grouped = king_names_basins.groupby('basin')
store = {}
for basin, glaciers in grouped:
	q = king_annual_glacier_flux.filter(items=glaciers.king_name).sum(axis=1)
	store[basin] = q
basin_q_king = pd.DataFrame(store)


## Import Rignot's by-basin measurements 
rignot_raw = pd.read_csv(folder_path + 'ice_discharge/Bamber+co-ords_updatedbasins.csv', 
	index_col='rignot_name')
rignot = rignot_raw.drop(labels=['x', 'y', 'area', 'F2000', 'F1996'], axis=1)
# Drop any entries without a basin
rignot = rignot[rignot.basin.notnull()]
rignot.basin = ['basin{:.0f}'.format(n) for n in rignot.basin]
# Sum up the few cases where there a multiple outlets per basin (e.g. Petermann)
rignot = rignot.groupby(rignot.basin).sum()
# Basins now unique and provide the row index
# Create a date index
date_ix = [pd.datetime(int(y), 1, 1) for y in rignot.columns]
# Transpose rows<-->columns
basin_q_rignot = rignot.T
# Rows are now dates, columns are basins
basin_q_rignot.index = date_ix


## Use King data to work out average % contribution of each glacier to basin's 
## outflow for pre-2000 segmentation
grouped = king_names_basins.groupby('basin')
store = {}
for basin, glaciers in grouped:	
	# Calculate mean annual basin-wide Q
	#  !! This is a duplicate of some code above
	q_basin = king_annual_glacier_flux.filter(items=glaciers.king_name).sum(axis=1).mean()
	# Calculate mean annual Q per glacier
	q_glaciers = king_annual_glacier_flux.filter(items=glaciers.king_name).mean(axis=0)
	# Calculate %contribution of each glacier to basin Q
	qp = (100 / q_basin) * q_glaciers
	# Store it
	store[basin] = qp
glacier_basin_contrib = pd.DataFrame(store)
# Collapse to glacier-basinid-contribution.
glacier_basin_contrib = pd.concat({
	'basin': glacier_basin_contrib.apply(lambda row: row.first_valid_index(), axis=1), 
	'contrib': glacier_basin_contrib.apply(lambda row: row.dropna().iloc[0], axis=1)
	}, axis=1)


## Scale Rignot basin data to remove offset relative to King
# Use only the common temporal period
# Basins without King data get set to NaN...
basin_scaling = basin_q_rignot['2000':'2009'].mean() - basin_q_king['2000':'2009'].mean()
basin_q_rignot_sc = basin_q_rignot - basin_scaling
# ...so substitute these back in using mean basin scaling
# n.b. some null columns remain because Rignot dataset does not
# contain data for every single basin that they defined.
basin_q_rignot_sc[basin_q_rignot_sc.isnull()] = basin_q_rignot - basin_scaling.mean()


## facet plot of basin-by-basin comparisons
if plot_figs:
	plt.figure(figsize=(10,7))
	n = 1
	for b in basin_q_rignot.columns:
		plt.subplot(5, 7, n)
		plt.title(b)
		if b in basin_q_king:
			plt.plot(basin_q_rignot.index, basin_q_rignot[b], 'blue')
			plt.plot(basin_q_king.index, basin_q_king[b], 'red')
			plt.plot(basin_q_rignot_sc.index, basin_q_rignot_sc[b], '--b')
			plt.xlim('2000-01-01', '2015-12-31')
			n += 1
	plt.tight_layout()
	plt.savefig(folder_path + 'outputs_Nov2017/figures/basins_comparison.pdf')	


## Split Rignot per-basin data out to King-defined outlets using King % contributions
# Produces a 'per-glacier' time series like King's
# We are trying to get a dataframe with glaciers=columns, years=rows
store = []
for basin in basin_q_rignot_sc.columns:

	# Time series of basin discharge
	basin_q = basin_q_rignot_sc.loc[:, basin]

	# Check to see if King has measurement in this basin
	#if basin in glacier_basin_contrib.columns:
	if basin in glacier_basin_contrib.basin.tolist():
		# Extract % contributions
		#contribs = glacier_basin_contrib.T.loc[basin, :]
		contribs = glacier_basin_contrib[glacier_basin_contrib.basin == basin]
		# Remove all glaciers which do not contribute
		contribs = contribs[contribs.notnull()]
		print(contribs)
		# Allocate year-by-year
		glacier_q = []
		for year, y_q in basin_q.iteritems():
			# Contribution of each glacier for one year
			g_y_q = (y_q / 100) *  contribs.contrib
			g_y_q.name = year
			glacier_q.append(g_y_q)
		# Dataframe with cols=glaciers, rows=years
		glaciers_ts = pd.DataFrame(glacier_q)	
		store.append(glaciers_ts)	

	# Retain basin as an outflow in its own right
	else:
		store.append(basin_q)

glaciers_rignot_sc = pd.concat(store, axis=1)
glaciers_rignot_sc = glaciers_rignot_sc.dropna(axis=1, how='all')


## Replace Rignot values with King or Enderlin from 2000 onward.
from copy import deepcopy
glaciers_combined = deepcopy(glaciers_rignot_sc)
glaciers_combined = glaciers_combined.reindex(pd.date_range('1958-01-01', '2016-01-01', freq='AS'))
for glacier in glaciers_rignot_sc.columns:
	# Insert King data where possible, 2000:2015
	if glacier in king.columns: 
		glaciers_combined.loc['2000':'2016', glacier] = king_annual_glacier_flux[glacier]
	else:
		print('NOTADDED: ' + glacier)


## Enderlin measured a few extra glaciers (basins) not in Rignot - add them in
# E.g. 'basin51' near Thule is identified in Rignot shapefile, but Rignot does
# not provide any discharge estimates for the basin.
# At time of this comment, all these only_enderlin glaciers are in basins 51, 52.
only_king = []
for glacier in king_annual_glacier_flux.columns:
	if glacier not in glaciers_combined.columns:
		glaciers_combined = pd.concat((glaciers_combined, king_annual_glacier_flux[glacier]), axis=1)
		only_king.append(king_annual_glacier_flux[glacier])
# Need to double-check for accidental duplicates		


## Optional visualisation of results so far...
# Load Jonathan's 2012 paper values for comparison
vals2012 = pd.read_csv(folder_path + 'Annual_fluxes_2012paper.txt',
	names=['Runoff', 'Tundra', 'Discharge', 'Total'], delim_whitespace=True)
vals2012.index = pd.date_range('1958-01-01', '2010-01-01', freq='1AS')

# Get only the rignot 'glaciers' which King has data for
rignot_only_king = glaciers_rignot_sc.filter(items=king.columns)

if plot_figs:
	plt.figure(figsize=(8,6))
	plt.plot(glaciers_combined.index, glaciers_combined.sum(axis=1), linewidth=4, marker='s', label='Combined (Rignot.Sc, Enderlin, King)', alpha=0.6)
	plt.plot(basin_q_rignot.index, basin_q_rignot.sum(axis=1), marker='x', label='Rignot', alpha=0.6)
	plt.plot(glaciers_rignot_sc.index, glaciers_rignot_sc.sum(axis=1), marker='+', label='Rignot scaled (all glaciers)', alpha=0.6)
	plt.plot(king_annual_glacier_flux.index, king_annual_glacier_flux.sum(axis=1), marker='^', label='king', alpha=0.6)
	plt.plot(rignot_only_king.index, rignot_only_king.sum(axis=1), marker='x', label='Rignot scaled (only for king glaciers)', alpha=0.6)
	plt.plot(vals2012.index, vals2012.Discharge, marker='*', label='Bamber2012', alpha=0.6)
	plt.legend()
	plt.savefig(folder_path + 'outputs_Nov2017/figures/discharge_ds_comparison.pdf')


## Correlation with ice-sheet-wide runoff

# Load routed runoff
runoff = xr.open_dataset(folder_path + 'outputs_Nov2017/FWF17_runoff_RACMO2.3p2.nc',
	chunks={'TIME':6})

# Load mask (as we're only interested in Greenland for this analysis)
Gr_land_filled = georaster.SingleBandRaster(folder_path + 'outputs_Nov2017/mask_LSMGr_filled.tif')

# Use routed runoff
annual_runoff = runoff.runoff_ice \
	.where(Gr_land_filled.r) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X', 'Y')) \
	.to_pandas()

# Running 5y average (y and previous 4 y)
# Set min_periods=1 so that 1958 retains a value, as we also have solid ice Q this year
runoff_5y = annual_runoff.rolling(5, min_periods=1).mean()

# Combine, and remove zero discharge values
runoff_discharge = pd.DataFrame({'runoff':runoff_5y, 'discharge':glaciers_combined.sum(axis=1)})
runoff_discharge[runoff_discharge.discharge == 0] = np.nan
runoff_discharge = runoff_discharge.dropna()
# don't use king values as they don't capture all discharge
runoff_discharge = runoff_discharge[:'2012']

runoff_discharge.to_csv(folder_path + 'outputs_Nov2017/runoff_discharge_values.csv')

# Define and fit model
X = runoff_discharge.runoff
y = runoff_discharge.discharge
X = sm.add_constant(X)
model = sm.OLS(y, X)
results = model.fit()
print(results.summary())
with open(folder_path + 'outputs_Nov2017/runoff_discharge_model.txt', 'w') as fh:
	print(results.summary(), file=fh)

if plot_figs:
	plt.figure(figsize=(6,4))
	runoff_discharge.plot(kind='scatter', x='runoff', y='discharge', marker='x', color='k')
	plt.plot(runoff_discharge.runoff, results.fittedvalues, '-b')
	plt.savefig(folder_path + 'outputs_Nov2017/figures/runoff_v_discharge.pdf')



# Errors for correlation / sigma values??

# Estimate whole time series solid ice discharge back in time
sid_est = results.predict(exog=sm.add_constant(runoff_5y))
sid = deepcopy(sid_est)
sid[runoff_discharge.index] = runoff_discharge.discharge

if plot_figs:
	plt.figure(figsize=(5,5))
	sid_est.plot(label='estimated')
	sid.plot(label='estimated+observed')
	sid.rolling(5, min_periods=1).mean().plot(label='e+o 5y')
	vals2012.Discharge.rolling(5, min_periods=1).mean().plot(label='Bamber2012')
	plt.legend()
	plt.ylim(160, 1250)
	plt.savefig(folder_path + 'outputs_Nov2017/figures/discharge_observed_v_estimated.pdf')


## Attribute ice-sheet-wide solid ice discharge to specific glaciers, using monthly contrib

# Rather than going via basin, we want to split into individual glaciers 
# directly, so we have some trickery to do...
mean_q = glaciers_combined.loc['2000':'2012'].mean()
perc_q = (1. / mean_q.sum()) * mean_q
# Distribute annual flux over all glaciers
sid_glaciers = sid.apply(lambda row: row * perc_q)
sid_glaciers[glaciers_combined.notnull()] = glaciers_combined

# We now have an annual-resolution time series for individual glaciers 1958-2012

## Now convert to monthly time series...

# First forward fill the annual values onto monthly
sid_glaciers_monthly = sid_glaciers.resample('1MS').ffill()
# Currently ends on 2016-01-01, push out to 2016-12-01 
sid_glaciers_monthly = sid_glaciers_monthly \
	.reindex(pd.date_range('1958-01-01', '2016-12-01', freq='1MS')) \
	.fillna(method='ffill')

# For non-king glaciers (the few remaining Rignot basins), now calculate average %contrib a month from King series
monthly_dist = monthly_perc.groupby(monthly_perc.index.month).mean().mean(axis=1) / 100
# Apply this scaling to the annual values
sid_glaciers_monthly_generic = sid_glaciers_monthly \
	.drop(labels=monthly_flux.columns, axis=1) \
	.apply(lambda row: monthly_dist.loc[row.name.month] * row, axis=1)

# Step 2: apply monthly glacier-specific distribution (pre-2000) or flux (2000+)
glac_monthly = monthly_perc.groupby(monthly_perc.index.month).mean() / 100
sid_glaciers_monthly_king = {}
for glacier in monthly_flux.columns:
	store = []
	for ix, val in sid_glaciers_monthly[glacier].iteritems():
		if ix < pd.datetime(2000, 1, 1):
			# Apply the mean monthly distribution to the runoff-estimated value
			store.append(glac_monthly[glacier].loc[ix.month] * val)
		else:
			# Insert the real data from the King dataset
			store.append(monthly_flux.loc[ix, glacier])
	sid_glaciers_monthly_king[glacier] = pd.Series(store, index=sid_glaciers_monthly_generic.index, name=glacier)

# Now concatenate King and non-King glaciers back into a single dataframe
sid_glaciers_monthly = pd.concat((sid_glaciers_monthly_generic, pd.DataFrame(sid_glaciers_monthly_king)), axis=1)

# Calculate final annual dataset
sid_glaciers_annual = sid_glaciers_monthly.resample('1AS').sum()




# Define efflux points for remaining basins
# These were picked manually using combined.shp in QGIS for guidance
store = {}
store['basin18'] = [-266731.308, -2711727.904]
store['basin2'] = [-84384.416, -908161.157]
store['basin3'] = [-4775.744, -911388.535]
store['basin35'] = [338940.078, -972708.729]
store['basin4'] = [279771.470, -881804.231]
xtra_basin_locs = pd.DataFrame(store).T
xtra_basin_locs.columns = ['x', 'y']
geometry = [Point(xy) for xy in zip(xtra_basin_locs.x, xtra_basin_locs.y)]
xtra_basins = gpd.GeoDataFrame({'king_name':xtra_basin_locs.index}, geometry=geometry, crs=pstere)


## Prepare coordinates for export
complete_outlet_points = pd.concat([king_pstere, xtra_basins], axis=0, ignore_index=True)
complete_outlet_points.index = complete_outlet_points.king_name
complete_outlet_points_geo = complete_outlet_points.to_crs({'init':'epsg:4326'})
# First sort alphabetically
complete_outlet_points_geo = complete_outlet_points_geo.sort_index(axis=0)
row_x = [r.x for r in complete_outlet_points_geo.geometry]
row_y = [r.y for r in complete_outlet_points_geo.geometry]
coords_df = pd.DataFrame(np.array([np.array(row_x), np.array(row_y)]), 
	index=['longitude_degrees', 'latitude_degrees'], 
	columns=complete_outlet_points_geo.index.values)

## Export CSV (glacier-by-glacier)
# First sort alphabetically
sid_glaciers_monthly = sid_glaciers_monthly.reindex_axis(sorted(sid_glaciers_monthly.columns), axis=1)
# Add lat/lon rows to top of frame
to_export = pd.concat((coords_df, sid_glaciers_monthly), axis=0)
# Export, to 1dp inline with source datasets
to_export.to_csv(folder_path + 'outputs_Nov2017/sid_glaciers_monthly_coords.csv',
	float_format='%.3f', date_format='%Y-%m-%d')



### ==========================================================================
### Gridding
grid = False

if grid:
	## The initial distance lookups are copied from arctic_runoff_routing.py
	print('Gridding . . . ')
	## Create look-up tree of coastline pixels
	# Load the 'distance from the ocean' raster
	dist_land = georaster.SingleBandRaster(folder_path + 'outputs_Nov2017/ROUTING_distances_landandice.tif')
	coast = np.where(dist_land.r == 1, 1, 0)

	# In pixel space, create a 2-d array of coastline pixels
	x = np.arange(0, dist_land.nx)
	y = np.arange(0, dist_land.ny)
	xi, yi = np.meshgrid(x, y)
	xp = xi[np.where(coast == 1, True, False)].flatten()
	yp = yi[np.where(coast == 1, True, False)].flatten()
	coast_points = np.zeros((len(xp), 2))
	coast_points[:, 0] = yp
	coast_points[:, 1] = xp

	# Create the coast lookup tree
	tree = spatial.cKDTree(coast_points)


	sid_grid = np.zeros((len(sid_glaciers_monthly), dist_land.ny, dist_land.nx))

	# Use 1dp rounded version of dataset
	sid_glaciers_monthly_round = round(sid_glaciers_monthly, 3)

	print('Iterating pixels . . . ')
	# Add SID from each location in turn
	for ix, row in complete_outlet_points.iterrows():
		# First convert to pixel coordinates
		x_px, y_px = dist_land.coord_to_px(row['geometry'].x, row['geometry'].y)
		# Now find nearest coastal pixel
		distance, index = tree.query((y_px, x_px), k=1)
		cpy, cpx = coast_points[index, :]
		# Some outlets discharge within same 5 km pixel, hence +=
		sid_grid[:, int(cpy), int(cpx)] += sid_glaciers_monthly_round[row['king_name']]

	print('Saving grids . . . ')
	## Convert to netCDF
	coords = {'TIME':runoff.TIME, 'Y':runoff.Y, 'X':runoff.X}
	# Convert to DataArray, integer
	sid = xr.DataArray(sid_grid, coords=coords, dims=['TIME', 'Y', 'X'],
		encoding={'dtype':'int16', 'scale_factor':0.001, 'zlib':True, '_FillValue':-9999})
	sid.name = 'Solid ice discharge'
	sid.attrs['long_name'] = 'Solid ice discharge'
	sid.attrs['units'] = 'km3'
	sid.attrs['grid_mapping'] = 'polar_stereographic'

	ds = xr.Dataset({'solid_ice':sid, 'lon':runoff.lon, 'lat':runoff.lat, 
		'polar_stereographic':runoff.polar_stereographic})

	# Main metadata
	ds.attrs['Conventions'] = 'CF-1.4'
	ds.attrs['history'] = 'This NetCDF generated using bitbucket atedstone/fwf/solid_ice_discharge.py using data provided by Michalea King, Ian Howat and Eric Rignot'
	ds.attrs['institution'] = 'University of Bristol (Andrew Tedstone)'
	ds.attrs['title'] = 'Monthly solid ice discharge from Greenland on a projected grid'

	# Additional geo-referencing
	ds.attrs['nx'] = float(dist_land.nx)
	ds.attrs['ny'] = float(dist_land.ny)
	ds.attrs['xmin'] = float(np.round(np.min(runoff.X), 0))
	ds.attrs['ymax'] = float(np.round(np.max(runoff.Y), 0))
	ds.attrs['spacing'] = 5000.

	# NC conventions metadata for dimensions variables
	ds.X.attrs['units'] = 'meters'
	ds.X.attrs['standard_name'] = 'projection_x_coordinate'
	ds.X.attrs['point_spacing'] = 'even'
	ds.X.attrs['axis'] = 'X'

	ds.Y.attrs['units'] = 'meters'
	ds.Y.attrs['standard_name'] = 'projection_y_coordinate'
	ds.Y.attrs['point_spacing'] = 'even'
	ds.Y.attrs['axis'] = 'Y'

	ds.TIME.attrs['standard_name'] = 'time'
	ds.TIME.attrs['axis'] = 'TIME'

	ds.to_netcdf(folder_path + 'outputs_Nov2017/FWF17_solidice_RACMO2.3p2.nc', format='NetCDF4')
