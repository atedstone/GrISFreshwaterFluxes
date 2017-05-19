"""
Generates a monthly solid ice discharge time series for the Greenland Ice 
Sheet, 1958 to 2015.

Uses three datasets: Rignot, Enderlin, King. 
- Rignot basins are scaled by Enderlin fluxes and used to provide data for 1958, 1964, 1992-1999 period.
- Enderlin is 2000-2012 'gold standard' dataset by-glacier, annual
- King is used to provide real monthly values where available, and otherwise to split annual fluxes into monthly by % distribution.

Andrew Tedstone, May 2017.
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import statsmodels.api as sm

import georaster

pstere = {'init':'epsg:3413'}
plot_figs = True

### Load ice discharge datasets

## King/Howat/OSU monthly dataset
# First create/tidy up the column headings
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

# Now use columns to load dataset
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
king_geo.to_file('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/king_locs.shp')
king_pstere = king_geo.to_crs(pstere)


## Enderlin annual dataset
# Specifying nrows cuts the various average rows off the bottom of the dataset 
enderlin_raw = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/GrIS_D-timeseries.txt', encoding='mac_roman', 
	delim_whitespace=True, index_col='Glacier_Name', nrows=178)

enderlin = enderlin_raw.drop(enderlin_raw.columns[0], axis=1)
enderlin = enderlin.drop(enderlin_raw.columns[1], axis=1)
enderlin = enderlin.T
enderlin.index = pd.date_range('2000-01-01', '2012-01-01', freq='1AS')
enderlin_cols = [c + '_discharge' for c in enderlin.columns]
enderlin.columns = enderlin_cols

# Enderlin glacier locations
enderlin_locs = pd.concat([enderlin_raw[enderlin_raw.columns[0]], enderlin_raw[enderlin_raw.columns[1]]], axis=1)
enderlin_locs = enderlin_locs.apply(pd.to_numeric, errors='ignore')
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

"""
 At this point I exported the king_equiv values above and then did some further
 manual re-mapping of a few labels.
 ...hence the manually-specified list of names below.
"""

 # List is same length as King dataset.
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
# Retain only discharge, get rid of std
king = king.filter(like='_discharge')
# Now remap King columns to Enderlin names
king.columns = enderlin_names_for_king

## Convert King's monthly rate to monthly flux then calculate total ice-sheet-wide annual flux
king = king.assign(year_length=np.where(king.index.is_leap_year, 366, 365))
king = king.assign(month_length=king.index.daysinmonth)
monthly_rate = king.filter(like='discharge')
monthly_flux = monthly_rate.apply(lambda x: (x / king.year_length.values) * king.month_length.values)
# Annual flux ice-sheet-wide
king_annual_total_flux = monthly_flux.resample('1AS').sum().T.sum()
# Annual flux per flacier
king_annual_glacier_flux = monthly_flux.resample('1AS').sum()

# Compare Enderlin and King
# This doesn't necessarily work currently
if plot_figs:
	plt.figure()
	n = 1
	for nk, ne in zip(king.filter(like='discharge').columns, enderlin_names_for_king):
		print(nk, ne)	
		plt.subplot(8, 6, n)
		plt.title(ne)
		plt.plot(enderlin.index, enderlin[ne], 'r')
		plt.plot(king_annual_glacier_flux.index, king_annual_glacier_flux[nk], 'b')
		n += 1


# Discharge record from two glaciers does not match - ukassorssuaq and torsukatat
# Nevertheless, continue for the moment
# Enderlin total for King glaciers:
enderlin_annual_subset = enderlin.filter(items=king.columns).T.sum()
print(enderlin_annual_subset)
print(king_annual_total_flux)


## Calculate % contribution of each month to annual flux...
# percentages make glaciers comparable to one another
# first convert to monthly perc of annual each year at each glacier - i.e. still a time series
monthly_perc = monthly_flux.groupby(monthly_flux.index.year).apply(lambda x: (100 / x.sum()) * x)
#monthly_perc.groupby(monthly_perc.index.month).mean()
## the monthly percentage is actually very similar every month on average, i.e. c. 8%/year.



## Load Rignot's drainage basins
rignot_basins = gpd.read_file('/home/at15963/Dropbox/work/papers/bamber_fwf/JLB_Analysis_2015/combined.shp')
# remove the 'basin' prefix
basins_pstere = rignot_basins.to_crs(pstere)
basins_pstere.columns = ['GRIDCODE', 'area', 'geometry', 'basin']


## Resolve Ellyns' individual outlet glaciers into their Rignot basins
# First do a spatial join, based on glacier point *within* basin poly
enderlin_names_basins = gpd.sjoin(enderlin_pstere, basins_pstere, how='left', op='within')

# Next, for the outlets that are not within a basin, resolve them to the nearest one
enderlin_names_basins = enderlin_names_basins.assign(dist=0)
unalloc = enderlin_names_basins[enderlin_names_basins.basin.isnull()]
for ix, row in unalloc.iterrows():
	dists = basins_pstere.distance(row.geometry)
	ix_min = dists.argmin()
	nearest_basin = basins_pstere.loc[ix_min, 'basin']
	enderlin_names_basins.loc[ix, 'basin'] = nearest_basin	
	enderlin_names_basins.loc[ix, 'dist'] = dists.min()


## For mass continuity with Rignot, move some sub-basins into bigger basins...
# allocate basins below Helheim into Helheim basin
enderlin_names_basins.loc[enderlin_names_basins.basin.str.contains('basin46'), 'basin'] = 'basin11'
enderlin_names_basins.loc[enderlin_names_basins.basin.str.contains('basin47'), 'basin'] = 'basin11'
# ikertivaq north, pamiataq, unnamed into ikertivaq (as sub-basin is on flow divide)
enderlin_names_basins.loc[enderlin_names_basins.basin.str.contains('basin48'), 'basin'] = 'basin12'
enderlin_names_basins.loc[enderlin_names_basins.basin.str.contains('basin49'), 'basin'] = 'basin12'
# Allocate all glaciers in basin38_* basins to basin38
enderlin_names_basins.loc[enderlin_names_basins.basin.str.contains('basin38_'), 'basin'] = 'basin38'


## Aggregate Ellyn's measurements to basin outputs
grouped = enderlin_names_basins.groupby('basin')
store = {}
for basin, glaciers in grouped:
	q = enderlin.filter(items=glaciers.enderlin_name).sum(axis=1)
	store[basin] = q
basin_q_enderlin = pd.DataFrame(store)


## Import Rignot's by-basin measurements 
rignot_raw = pd.read_excel('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/Bamber+co-ords_updatedbasins.xlsx', 
	index_col='rignot_name')
rignot = rignot_raw.drop(labels=['x', 'y', 'area', 'F2000', 'F1996'], axis=1)
# Drop any entries without a basin
rignot = rignot[rignot.basin.notnull()]
rignot.basin = ['basin{:.0f}'.format(n) for n in rignot.basin]
# Sum up the few cases where there a multiple outlets per basin (e.g. Petermann)
rignot = rignot.groupby(rignot.basin).sum()
# Basins now unique and provide the row index
# Create a date index
date_ix = [pd.datetime(y, 1, 1) for y in rignot.columns]
# Transpose rows<-->columns
basin_q_rignot = rignot.T
# Rows are now dates, columns are basins
basin_q_rignot.index = date_ix


## Use Enderlin data to work out % contribution of each glacier to basin's outflow
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
# Collapse to glacier-basinid-contribution.
#basin_ids = glacier_basin_contrib.apply(lambda row: row.first_valid_index(), axis=1)
#basin_contrib = glacier_basin_contrib.apply(lambda row: row.dropna().iloc[0], axis=1)
glacier_basin_contrib = pd.concat({
	'basin': glacier_basin_contrib.apply(lambda row: row.first_valid_index(), axis=1), 
	'contrib': glacier_basin_contrib.apply(lambda row: row.dropna().iloc[0], axis=1)
	}, axis=1)


## Scale Rignot basin data to remove offset relative to Enderlin
# Use only the common temporal period
# Basins without Enderlin data get set to NaN...
basin_scaling = basin_q_rignot['2000':'2009'].mean() - basin_q_enderlin['2000':'2009'].mean()
basin_q_rignot_sc = basin_q_rignot - basin_scaling
# ...so substitute these back in using mean basin scaling
# n.b. some null columns remain because Rignot dataset does not
# contain data for every single basin that they defined.
basin_q_rignot_sc[basin_q_rignot_sc.isnull()] = basin_q_rignot - basin_scaling.mean()


## facet plot of basin-by-basin comparisons
if plot_figs:
	plt.figure()
	n = 1
	for b in basin_q_rignot.columns:
		plt.subplot(8, 6, n)
		plt.title(b)
		try:
			plt.plot(basin_q_rignot.index, basin_q_rignot[b], 'blue')
			plt.plot(basin_q_enderlin.index, basin_q_enderlin[b], 'red')
			plt.plot(basin_q_rignot_sc.index, basin_q_rignot_sc[b], '--b')
			plt.xlim('2000-01-01', '2015-12-31')
		except KeyError:
			# This logic means we don't see basins that only have Enderlin data
			pass
		n += 1


## Split Rignot per-basin data out to Enderlin-defined outlets using Enderlin % contributions
# Produces a 'per-glacier' time series like Enderlin's
# We are trying to get a dataframe with glaciers=columns, years=rows
store = []
for basin in basin_q_rignot_sc.columns:

	# Time series of basin discharge
	basin_q = basin_q_rignot_sc.loc[:, basin]

	# Check to see if Enderlin has measurement in this basin
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
glaciers_combined = glaciers_combined.reindex(pd.date_range('1958-01-01', '2015-01-01', freq='AS'))
for glacier in glaciers_rignot_sc.columns:
	# Insert King data where possible, 2000:2015
	if glacier in king.columns: ## Not using the correct names at the moment!
		glaciers_combined.loc['2000':'2015', glacier] = king_annual_glacier_flux[glacier]
		pass
	# Otherwise insert Enderlin data where it exists
	elif glacier in enderlin.columns:
		glaciers_combined.loc['2000':'2012', glacier] = enderlin[glacier]


## Enderlin measured a few extra glaciers (basins) not in Rignot - add them in
# E.g. 'basin51' near Thule is identified in Rignot shapefile, but Rignot does
# not provide any discharge estimates for the basin.
# At time of this comment, all these only_enderlin glaciers are in basins 51, 52.
only_enderlin = []
for glacier in enderlin.columns:
	if glacier not in glaciers_combined.columns:
		glaciers_combined = pd.concat((glaciers_combined, enderlin[glacier]), axis=1)
		only_enderlin.append(enderlin[glacier])
# Need to double-check for accidental duplicates		


## Optional visualisation of results so far...
# Load Jonathan's 2012 paper values for comparison
vals2012 = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/Annual_fluxes_2012paper.txt',
	names=['Runoff', 'Tundra', 'Discharge', 'Total'], delim_whitespace=True)
vals2012.index = pd.date_range('1958-01-01', '2010-01-01', freq='1AS')

# Get only the rignot 'glaciers' which Enderlin has data for
rignot_only_enderlin = glaciers_rignot_sc.filter(items=enderlin.columns)

if plot_figs:
	plt.figure()
	plt.plot(glaciers_combined.index, glaciers_combined.sum(axis=1), linewidth=4, marker='s', label='Combined (Rignot.Sc, Enderlin, King)', alpha=0.6)
	plt.plot(basin_q_rignot.index, basin_q_rignot.sum(axis=1), marker='x', label='Rignot', alpha=0.6)
	plt.plot(glaciers_rignot_sc.index, glaciers_rignot_sc.sum(axis=1), marker='+', label='Rignot scaled (all glaciers)', alpha=0.6)
	plt.plot(enderlin.index, enderlin.sum(axis=1), marker='^', label='Enderlin', alpha=0.6)
	plt.plot(rignot_only_enderlin.index, rignot_only_enderlin.sum(axis=1), marker='x', label='Rignot scaled (only for Enderlin glaciers)', alpha=0.6)
	plt.plot(vals2012.index, vals2012.Discharge, marker='*', label='Bamber2012', alpha=0.6)
	plt.legend()


## Correlation with ice-sheet-wide runoff

# Load runoff
runoff = xr.open_dataset('/scratch/process/RACMO2.3_GRN11_runoff_monthly_1958-2015_pstere.nc')

# Load mask (as we're only interested in Greenland for this analysis)
masks = xr.open_dataset('/scratch/process/RACMO2.3_GRN11_masks_pstere.nc')
GrIS_mask = masks.GrIS_mask

# Convert from mm w.e. to km3
runoff_flux = runoff.runoff * (5*5) / 1.0e6
annual_runoff = runoff_flux \
	.where(GrIS_mask) \
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

# Define and fit model
X = runoff_discharge.runoff
y = runoff_discharge.discharge
X = sm.add_constant(X)
model = sm.OLS(y, X)
results = model.fit()
print(results.summary())

if plot_figs:
	plt.figure()
	runoff_discharge.plot(kind='scatter', x='runoff', y='discharge', marker='x', color='k')
	plt.plot(runoff_discharge.runoff, results.fittedvalues, '-b')


# # try a rignot-only, a la Bamber 2012:
# # we wouldn't expect precise match as Rignot values have been scaled by Enderlin's...
# runoff_discharge_r = pd.DataFrame({'runoff':runoff_5y, 'discharge':basin_q_rignot.sum(axis=1)})
# runoff_discharge_r[runoff_discharge_r.discharge == 0] = np.nan
# runoff_discharge_r = runoff_discharge_r.dropna()
# X = runoff_discharge_r.runoff
# y = runoff_discharge_r.discharge
# X = sm.add_constant(X)
# model = sm.OLS(y, X)
# results = model.fit()
# print(results.summary())

# # check correlation of bamber d/s
# X = vals2012.Runoff.rolling(5, min_periods=1).mean()
# y = vals2012.Discharge
# X = sm.add_constant(X)
# model = sm.OLS(y, X)
# results = model.fit()
# print(results.summary())




# Errors for correlation / sigma values??

# Estimate whole time series solid ice discharge back in time
sid_est = results.predict(exog=sm.add_constant(runoff_5y))
sid = deepcopy(sid_est)
sid[runoff_discharge.index] = runoff_discharge.discharge

if plot_figs:
	plt.figure()
	sid_est.plot(label='estimated')
	sid.plot(label='estimated+observed')
	sid.rolling(5, min_periods=1).mean().plot(label='e+o 5y')
	vals2012.Discharge.rolling(5, min_periods=1).mean().plot(label='Bamber2012')
	plt.legend()
	plt.ylim(160, 1250)


## Attribute ice-sheet-wide solid ice discharge to specific glaciers, using monthly contrib

# Rather than going via basin, we want to split into individual glaciers 
# directly, so we have some trickery to do...
mean_q = glaciers_combined.loc['2000':'2012'].mean()
perc_q = (1. / mean_q.sum()) * mean_q
# Distribute annual flux over all glaciers
sid_glaciers = sid.apply(lambda row: row * perc_q)
sid_glaciers[glaciers_combined.notnull()] = glaciers_combined

# We now have an annual-resolution time series for individual glaciers 1958-2012
# Values for 2013-2015 are not yet correct.


## Estimate solid ice discharge forward in time (up to 2015)
"""
First deal with non-King outlets. Drop them from percentage contribs look up
and then use remaining glaciers to allocate flux.
Then essentially allocate the rest of the flux by adding King outlets back on.

The underlying rationale here is that we know the %contrib which each King 
outlet makes based on our analysis above. This approach doesn't preserve the 
mass of the modelled ice-sheet-wide SID for 2013 to 2015 but DOES preserve the
mass output by the King outlets.
"""
perc_q_remaining = perc_q.drop(labels=king_annual_glacier_flux.columns)
sid_1315 = sid['2013':'2015'].apply(lambda row: row * perc_q_remaining)
# Here we bash the King estimates on...not strictly needed as we fill these data in below, monthly.
sid_1315 = pd.concat((sid_1315, king_annual_glacier_flux['2013':'2015']), axis='columns')

sid_glaciers = sid_glaciers.drop(labels=pd.date_range('2013-01-01', '2015-01-01', freq='AS'), axis=0)
sid_glaciers = pd.concat((sid_glaciers, sid_1315))


## Now convert to monthly time series...

# First forward fill the annual values onto monthly
sid_glaciers_monthly = sid_glaciers.resample('1MS').ffill()
# Currently ends on 2015-01-01, push out to 2015-12-01 
sid_glaciers_monthly = sid_glaciers_monthly \
	.reindex(pd.date_range('1958-01-01', '2015-12-01', freq='1MS')) \
	.fillna(method='ffill')

# For non-king glaciers, now calculate average %contrib a month from King series
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

# Export
sid_glaciers_monthly.to_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/sid_glaciers_monthly.csv')


if plot_figs:
	plt.figure()
	plt.plot(glaciers_combined.index, glaciers_combined.sum(axis=1), linewidth=4, marker='s', label='Combined (Rignot.Sc, Enderlin, King)', alpha=0.6)
	plt.plot(basin_q_rignot.index, basin_q_rignot.sum(axis=1), marker='x', label='Rignot', alpha=0.6)
	plt.plot(glaciers_rignot_sc.index, glaciers_rignot_sc.sum(axis=1), marker='+', label='Rignot scaled (all glaciers)', alpha=0.6)
	plt.plot(enderlin.index, enderlin.sum(axis=1), marker='^', label='Enderlin', alpha=0.6)
	plt.plot(rignot_only_enderlin.index, rignot_only_enderlin.sum(axis=1), marker='x', label='Rignot scaled (only for Enderlin glaciers)', alpha=0.6)
	plt.plot(vals2012.index, vals2012.Discharge, marker='*', label='Bamber2012', alpha=0.6)
	plt.plot(sid_glaciers_annual.index, sid_glaciers_annual.sum(axis=1), marker='*', label='final monthly agg. to annual', alpha=0.6)
	plt.legend()


# Calculate an efflux point for basins (i.e. those without defined outlet glaciers)
basins_pstere.iloc[2].geometry.exterior.coords.xy



# Need to generate a complete coordinate series of all glacier outlet points

# What to do about basins that only have Enderlin values?
# Basins which only have Rignot are quite easy as Rignot goes back to 1992. (although what about forward?)



# Then can calculate ice-sheet wide Q 1992-2015 and produce correlation with runoff.


# Now deal with pre-1992 data set...
## Take a very similar approach with annual Q estimates from correlation
# but rather than thinking in basin terms, just attribute direct to each outlet
# some outlets also have their own seasonal cycle to consider (c.f. King)






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




















"""
##Previous attempts to join all three time series together using glacier names (NOT containing basins)

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
"""

"""
##old rignot basin merging logic


rignot_names_basins = pd.Series(rignot.index)
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
"""