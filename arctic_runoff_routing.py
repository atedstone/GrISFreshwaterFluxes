"""
arctic_runoff_routing.py

Do standard hydrological routing to identify points where ice runoff meets
ocean (land?? in case of wGrIS - to be resolved)

Route monthly runoff to outflux points


@author Andrew Tedstone a.j.tedstone@bristol.ac.uk
@created 2016-December-08

"""

import numpy as np
from pygeoprocessing import routing
from skimage import morphology
import xarray as xr
import pandas as pd

import georaster

dem_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/topo_racmo_EPSG3413_5km.tif'

## Step 0: fill pits
dempits_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/routing/racmo_EPSG3413_5km_0dempits.tif'
routing.fill_pits(dem_uri, dempits_uri)

## Step 1: flow direction
flowdir_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/routing/racmo_EPSG3413_5km_1flowdir.tif'
routing.flow_direction_d_inf(dempits_uri, flowdir_uri)

## Step 2: flow accumulation
accum_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/routing/racmo_EPSG3413_5km_2accum.tif'
routing.flow_accumulation(flowdir_uri, dempits_uri, accum_uri)


## Create an absorption raster, in this case absorption is zero
# Use the accum_uri for georeferncing
absorp_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/routing/racmo_EPSG3413_5km_absorp.tif'
absorp = georaster.SingleBandRaster(accum_uri)
zeros = np.zeros(absorp.r.shape)
absorp.r = zeros
absorp.save_geotiff(absorp_uri, dtype=georaster.gdal.GDT_UInt16)
absorp = None


## Identify outflux pixels (ICE MASK)
ice_mask_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/mask_racmo_EPSG3413_5km.tif'
ice_mask = georaster.SingleBandRaster(ice_mask_uri)
# Grab geo-referencing, for writing out new geoTIFFs
trans = ice_mask.trans
proj4 = ice_mask.srs.ExportToProj4()
# Get distances
out, distances = morphology.medial_axis(ice_mask.r, return_distance=True)
out_locations = (distances == 1)
georaster.simple_write_geotiff(
	'/media/sf_scratch/distances.tif',
	distances,
	trans, 
	proj4=proj4,
	dtype=georaster.gdal.GDT_UInt16
	)


## Step 3: route monthly runoff

runoff = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_EPSG3413_runoff_monthly_1958.nc')


store = np.zeros((12, ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
n = 0
for date in pd.date_range('1958-01-01', '1958-12-01', freq='1MS'):
#for date in runoff.TIME:
	print(date)
	r_month = runoff.runoff.sel(TIME=date)
	r_month = r_month.where(np.flipud(ice_mask.r) == 1)
	# Convert mmWE to flux per grid cell
	r_month = r_month * (5*5) / 1.0e6
	georaster.simple_write_geotiff(
		'/media/sf_scratch/TMP_runoff_for_month.tif',
		np.flipud(r_month.values),
		trans,
		proj4=proj4,
		dtype=georaster.gdal.GDT_Float64
		)
	# output a temporary geotiff, simplest way of supplying to routing??
	# r_month_mem = georaster.SingleBandRaster.from_array(
	# 		np.flipud(r_month.values),
	# 		trans, 
	# 		proj4,
	# 		gdal_dtype=georaster.gdal.GDT_Float64
	# 		)
	# r_month_mem.save_geotiff('/media/sf_scratch/TMP_runoff_for_month.tif', dtype=georaster.gdal.GDT_Float64)
	# r_month_mem = None

	# For runoff, could either reference direct into netCDF (potentially tricky)
	#source_uri = '/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/racmo_pstere_5km_runoff_1958-07.tif'
	source_uri = '/media/sf_scratch/TMP_runoff_for_month.tif'
	loss_uri = '/media/sf_scratch/fwf_routing/loss.tif'
	flux_uri = '/media/sf_scratch/fwf_routing/flux.tif'
	routing.route_flux(
			flowdir_uri, dempits_uri, source_uri, absorp_uri,
			loss_uri, flux_uri, 'flux_only')

	f = georaster.SingleBandRaster(flux_uri)
	store[n,:,:] = f.r
	n += 1
	# MUST close handle to flux file otherwise pygeoprocessing can't overwrite on next iteration
	f = None

