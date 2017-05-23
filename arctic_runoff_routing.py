"""
arctic_runoff_routing.py

Uses gridded RACMO runoff, masks and topogrphy created using reproject_racmo.py.

Routes runoff from ice sheets and caps to their margins, using the 
pygeoprocessing toolbox.

Water is moved from the margin to the coast along the euclidean distance, 
determined using a KD-tree.

@author Andrew Tedstone a.j.tedstone@bristol.ac.uk
@created 2016-December-08

"""

import numpy as np
from pygeoprocessing import routing
from skimage import morphology
import xarray as xr
import pandas as pd
from scipy import spatial 
import os
from osgeo import gdal

import georaster

# If true, only process 1958
debug = False

PROCESS_DIR = os.environ['PROCESS_DIR']

# ----------------------------------------------------------------------------

## Only use topography of ice surfaces - route to ice sheet + cap margins.
dem_uri = PROCESS_DIR + 'dem_ice_racmo_EPSG3413_5km.tif'
dem_sm_uri = PROCESS_DIR + 'dem_ice_sm_racmo_EPSG3413_5km.tif'
# Set nans to zero - I think this is needed for pit filling?
dem = georaster.SingleBandRaster(dem_uri)
dem.r = np.where(np.isnan(dem.r), 0, dem.r)
dem.save_geotiff(dem_sm_uri)
dem = None


## Step 0: fill pits
dempits_uri = PROCESS_DIR + 'ROUTING_racmo_EPSG3413_5km_dem_filledpits.tif'
routing.fill_pits(dem_sm_uri, dempits_uri)


## Step 1: flow direction
flowdir_uri = PROCESS_DIR + 'ROUTING_racmo_EPSG3413_5km_1flowdir.tif'
routing.flow_direction_d_inf(dempits_uri, flowdir_uri)


## Step 2: flow accumulation
accum_uri = PROCESS_DIR + 'ROUTING_racmo_EPSG3413_5km_2accum.tif'
routing.flow_accumulation(flowdir_uri, dempits_uri, accum_uri)


## Create an absorption raster, in this case absorption is zero
# We need an absorption raster for the final routing step to work
# Use the accum_uri for georeferncing
absorp_uri = PROCESS_DIR + 'ROUTING_racmo_EPSG3413_5km_absorp.tif'
absorp = georaster.SingleBandRaster(accum_uri)
zeros = np.zeros(absorp.r.shape)
absorp.r = zeros
absorp.save_geotiff(absorp_uri, dtype=gdal.GDT_UInt16)
absorp = None


## Create distance raster of ice mask in order to find ice sheet edges later
# We are also loading the ice mask to use later in order to mask out tundra runoff
ice_mask_uri = PROCESS_DIR + 'mask_ice_racmo_EPSG3413_5km.tif'
ice_mask = georaster.SingleBandRaster(ice_mask_uri)
# Grab geo-referencing, for writing out new geoTIFFs
trans = ice_mask.trans
proj4 = ice_mask.srs.ExportToProj4()
# Get distances
out, distances = morphology.medial_axis(ice_mask.r, return_distance=True)
georaster.simple_write_geotiff(
	PROCESS_DIR + 'ROUTING_distances_ice.tif',
	distances,
	trans, 
	proj4=proj4,
	dtype=gdal.GDT_UInt16
	)


## Create distance raster of land mask so we can find coast later.
landice_mask_uri = PROCESS_DIR + 'mask_landandice_racmo_EPSG3413_5km.tif'
landice_mask = georaster.SingleBandRaster(landice_mask_uri)
# Grab geo-referencing, for writing out new geoTIFFs
trans = landice_mask.trans
proj4 = landice_mask.srs.ExportToProj4()
# Get distances
out, distances = morphology.medial_axis(landice_mask.r, return_distance=True)
georaster.simple_write_geotiff(
	PROCESS_DIR + 'ROUTING_distances_landandice.tif',
	distances,
	trans, 
	proj4=proj4,
	dtype=gdal.GDT_UInt16
	)



### Logic to move water from ice sheet margins to the coast cells

## Create look-up tree of coastline pixels
# Load the 'distance from the ocean' raster
dist_land = georaster.SingleBandRaster(PROCESS_DIR + 'ROUTING_distances_landandice.tif')
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

## Find the ice sheet margins
# Essentially 'Distance from land'
ice_dist = georaster.SingleBandRaster(PROCESS_DIR + 'ROUTING_distances_ice.tif')
ice = np.where(ice_dist.r == 1, 1, 0)

## Create DEM with ice margin elevations set to zero
# for mass conservation
dempits = georaster.SingleBandRaster(dempits_uri)
ice_dem_zeros = np.where(ice, 0, dempits.r)
georaster.simple_write_geotiff(
	PROCESS_DIR + 'ROUTING_dempits_zeromargins.tif',
	ice_dem_zeros,
	trans, 
	proj4=proj4,
	dtype=gdal.GDT_UInt16
	)

## In pixel space, create a 2-d array of ice margin pixels
x = np.arange(0, ice_dist.nx)
y = np.arange(0, ice_dist.ny)
xi, yi = np.meshgrid(x, y)
xp = xi[ice == 1].flatten()
yp = yi[ice == 1].flatten()
ice_points = np.zeros((len(xp), 2))
ice_points[:, 0] = yp
ice_points[:, 1] = xp

## Load land mask to route tundra runoff (euclidean distance)
land_mask_uri = PROCESS_DIR + 'mask_land_racmo_EPSG3413_5km.tif'
land_mask = georaster.SingleBandRaster(land_mask_uri)
x = np.arange(0, land_mask.nx)
y = np.arange(0, land_mask.ny)
xi, yi = np.meshgrid(x, y)
xp = xi[land_mask.r == 1].flatten()
yp = yi[land_mask.r == 1].flatten()
land_points = np.zeros((len(xp), 2))
land_points[:, 0] = yp
land_points[:, 1] = xp

## Step 3: route monthly runoff
# Open the gridded runoff file          
runoff = xr.open_dataset(PROCESS_DIR + 'RACMO2.3_GRN11_runoff_monthly_1958-2015_pstere.nc')

if debug == True:
	store_ice = np.zeros((12*7, ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	store_tundra = np.zeros((12*7, ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	dates = pd.date_range('1958-01-01', '1964-12-01', freq='1MS')
else:
	dates = runoff.TIME.to_index()
	store_ice = np.zeros((len(dates), ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	store_tundra = np.zeros((len(dates), ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))

n = 0
ice_r_pre = []
ice_r_post = []
for date in dates:
	print(date)
	# Import runoff data grid from netcdf
	r_month = np.flipud(runoff.runoff.sel(TIME=date).values.squeeze())
	# Convert mmWE to flux per grid cell
	r_month = np.where(r_month > 0, r_month * (5*5) / 1.0e6, 0)
	# Get ice component of runoff
	r_month_ice = np.where(ice_mask.r == 1, r_month, np.nan)
	ice_r_pre.append(np.nansum(r_month_ice))
	# Save to a geoTIFF so that the routing toolbox can ingest it
	georaster.simple_write_geotiff(
		PROCESS_DIR + 'TMP_runoff_for_month_%s.tif' % n,
		r_month_ice,
		trans,
		proj4=proj4,
		dtype=gdal.GDT_Float32
		)

	# Route the runoff to the ice margins
	source_uri = PROCESS_DIR + 'TMP_runoff_for_month_%s.tif' % n
	loss_uri = PROCESS_DIR + 'loss.tif'
	flux_uri = PROCESS_DIR + 'flux_month%s.tif' % n
	routing.route_flux(
			flowdir_uri, PROCESS_DIR + 'ROUTING_dempits_zeromargins.tif', source_uri, absorp_uri,
			loss_uri, flux_uri, 'flux_only')

	# Open up the fluxes so we can route from ice margin to coast.
	flux = georaster.SingleBandRaster(flux_uri)

	# Assign the runoff from each ice margin pixel to its nearest coastal pixel.
	# Uses the look-up tree that we created above.
	# Create the grid to save the coastal outflux values to
	ice_r_post.append(np.nansum(flux.r[ice == 1]))
	coast_grid_ice = np.zeros((dist_land.ny, dist_land.nx))
	for ry, rx in ice_points:
		distance, index = tree.query((ry, rx), k=1)
		cpy, cpx = coast_points[index, :]
		coast_grid_ice[int(cpy), int(cpx)] += flux.r[int(ry), int(rx)]	

	# Save coastal fluxes to store
	store_ice[n,:,:] = np.flipud(coast_grid_ice)

	# MUST close handle to flux file otherwise pygeoprocessing can't overwrite on next iteration
	f = None


	## Route tundra fluxes too ...
	# Do by euclidean distance
	r_month_tundra = np.where(land_mask.r == 1, r_month, np.nan)
	coast_grid_tundra = np.zeros((dist_land.ny, dist_land.nx))
	for ry, rx in land_points:
		distance, index = tree.query((ry, rx), k=1)
		cpy, cpx = coast_points[index, :]
		coast_grid_tundra[int(cpy), int(cpx)] += r_month_tundra[int(ry), int(rx)]

	# Save coastal fluxes to store
	store_tundra[n,:,:] = np.flipud(coast_grid_tundra)

	n += 1


"""
Masks to include in this product:
* Greenland              >
* Canadian High Arctic   > Could produce single combined mask with numbers for each area
* Svalbard               >
* Iceland                >

* Ocean basins
"""

#routed_old = xr.open_dataset('~/Dropbox/RACMO23_routed_1958_2015.nc')
#store_ice = routed_old.runoff_ice.to_masked_array()
#store_tundra = routed_old.runoff_tundra.to_masked_array()

coords = {'TIME':runoff.TIME, 'Y':runoff.Y, 'X':runoff.X}

da_ice = xr.DataArray(np.round(store_ice, 2) * 100, 
	coords=coords, 
	dims=['TIME', 'Y', 'X'], 
	encoding={'dtype':'int16', 'scale_factor':0.01, 'zlib':True})
da_ice.name = 'Ice sheet runoff'
da_ice.attrs['long_name'] = 'Ice sheet runoff'
da_ice.attrs['units'] = 'km3'
da_ice.attrs['grid_mapping'] = 'polar_stereographic'

da_tundra = xr.DataArray(np.round(store_tundra, 2) * 100, 
	coords=coords, 
	dims=['TIME', 'Y', 'X'], 
	encoding={'dtype':'int16', 'scale_factor':0.01, 'zlib':True})
da_tundra.name = 'Tundra runoff'
da_tundra.attrs['long_name'] = 'Tundra runoff'
da_tundra.attrs['units'] = 'km3'
da_tundra.attrs['grid_mapping'] = 'polar_stereographic'

ds = xr.Dataset({'runoff_ice':da_ice, 
	'runoff_tundra':da_tundra,
	'lon':runoff.lon,
	'lat':runoff.lat,
	'polar_stereographic':runoff.polar_stereographic})

# Main metadata
ds.attrs['Conventions'] = 'CF-1.4'
ds.attrs['history'] = 'This NetCDF generated using bitbucket atedstone/fwf/arctic_runoff_routing.py using data output by fwf/reproject_racmo.py\n' + runoff.history
ds.attrs['institution'] = 'University of Bristol (Andrew Tedstone)'
ds.attrs['title'] = 'Monthly ice sheet and tundra runoff routed to coastal pixels'

# Additional geo-referencing
ds.attrs['nx'] = float(dempits.nx)
ds.attrs['ny'] = float(dempits.ny)
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

ds.to_netcdf('/home/at15963/Dropbox/work/papers/bamber_fwf/outputs/FWF17_runoff.nc', format='NetCDF4')

pre = pd.Series(ice_r_pre, index=dates)
pre.to_csv(PROCESS_DIR + 'runoff_monthly_totals_pre_routing.csv')
post = pd.Series(ice_r_post, index=dates)
post.to_csv(PROCESS_DIR + 'runoff_monthly_totals_post_routing.csv')
