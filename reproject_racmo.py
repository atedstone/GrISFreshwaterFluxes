"""
Project RACMO 2.3 data onto a regular grid using interpolation.
Outputs an equivalent NetCDF of the time series.

Remember to set environment variable PROCESS_DIR, with trailing slash

example for polar stereographic full time series:
reproject_racmo.py /scratch/L0data/RACMO/RACMO2.3_GRN11_runoff_monthly_1958-2015.nc /scratch/L0data/RACMO/RACMO2.3_GRN11_masks.nc 1958-01-01 2015-12-31 

@author Andrew Tedstone a.j.tedstone@bristol.ac.uk
@created 2016-December-05

"""

import xarray as xr
import pyproj
from scipy import interpolate
import pandas as pd
import numpy as np
from osgeo import osr, gdal
from scipy import ndimage
import os
import argparse
import datetime as dt
import math
import georaster

parser = argparse.ArgumentParser(description='Project geographic RACMO data for Arctic domain')

# Positional arguments
parser.add_argument('fn_RACMO', type=str, help='str, absolute path to RACMO netCDF file')
parser.add_argument('fn_RACMO_masks', type=str, help='str, absolute path to RACMO masks netCDF file')
parser.add_argument('date_start', type=str, help='str, start date of RACMO data in format yyyy-mm-dd')
parser.add_argument('date_end', type=str, help='str, end date of RACMO data in format yyyy-mm-dd')

parser.add_argument('-grid', type=str, dest='grid', default='pstere', help='str, grid to use, pstere or bamber')
parser.add_argument('-test', dest='test', action='store_true', help='Test script (run only for first year of data)')

args = parser.parse_args()

##############################################################################

# Set grid to 'pstere' or 'bamber'
# 'bamber' grid is GrIS only, 'pstere' will result in grid across full RACMO extent
# Use 'bamber' grid only for testing to ensure that interpolation result matches 2012 FWF paper results
# Characteristics of each grid are hard-coded into logic below
grid = args.grid

PROCESS_DIR = os.environ['PROCESS_DIR']

# ============================================================================

def pstere_lat2k0(lat):
	""" Given lat in degrees, return k0 

	Implementation of eq. 1 in:
	Rollins, 2011. Computation of scale-factor and standard parallel for the 
	polar stereographic projection.
	URL: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
	Last retrieved 2012-Dec-07.

	"""

	ecc = 0.0818191908426215
	lat = np.deg2rad(70)
	k90 = np.sqrt(np.power((1. + ecc), (1. + ecc)) * np.power((1. - ecc), (1. - ecc)))
	k0_1 = ((1. + np.sin(lat)) / 2.) 
	k0_2 = (k90 / np.sqrt(np.power(1. + (ecc * np.sin(lat)), (1.+ecc)) * np.power(1. - (ecc * np.sin(lat)), (1.-ecc))))
	k0 = k0_1 * k0_2

	return k0


# ============================================================================


## Determine grid projection
# 5 km grid is hard-coded
if args.grid == 'pstere':
	grid_proj = pyproj.Proj('+init=EPSG:3413')
	fn_mask_ice = PROCESS_DIR + 'mask_ice_racmo_EPSG3413_5km.tif'
	fn_mask_land = PROCESS_DIR + 'mask_land_racmo_EPSG3413_5km.tif'
	fn_mask_greenland = PROCESS_DIR + 'mask_greenland_racmo_EPSG3413_5km.tif'
	fn_dem_ice = PROCESS_DIR + 'dem_ice_racmo_EPSG3413_5km.tif'
	fn_mask_landandice = PROCESS_DIR + 'mask_landandice_racmo_EPSG3413_5km.tif'
elif args.grid == 'bamber':
	grid_proj = pyproj.Proj('+proj=sterea +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

# Create SRS representation
srs = osr.SpatialReference()
srs.ImportFromProj4(grid_proj.srs)


## Open RACMO runoff dataset
# Times aren't encoded in format understandable by xarray
runoff = xr.open_dataset(args.fn_RACMO,
		decode_times=False, chunks={'time':15})
# Apply time dimension in format understandable by xarray
times = pd.date_range(args.date_start, args.date_end, freq='1MS')
runoff['time'] = times


## Extract and then project lat/lon coordinates to grid projection
lon = runoff.LON.values.flatten()
lat = runoff.LAT.values.flatten()
latlon = np.vstack((lat, lon)).T
(x, y) = grid_proj(lon, lat)
xy = np.vstack((x, y)).T

def grid_ceil(n):
	return math.ceil(n / 10000) * 10000

def grid_floor(n):
	return math.floor(n / 10000) * 10000

## Create grid at 5 km resolution
if grid == 'pstere':
	yi, xi = np.mgrid[grid_floor(np.min(y)):grid_ceil(np.max(y)):5000, grid_floor(np.min(x)):grid_ceil(np.max(x)):5000]
elif grid == 'bamber':
	yi, xi = np.mgrid[-3400000:-600000:5000, -800000:700000:5000]


## Calculate scaling factors for grid (after Bamber Interpolate_racmo.BAK)
# Calculate lat/lon coordinates for each point of projected grid
(grid_lon, grid_lat) = grid_proj(xi.flatten(), yi.flatten(), inverse=True)
# Convert to paired columns
grid_latlon = np.vstack((grid_lat, grid_lon)).T
# Convert to radians
grid_latlon_radians = grid_latlon / 57.29578
# Calculate latitudinal scaling
m = 2.0 * np.tan(45.0 / 57.29578 - (grid_latlon_radians[:, 0] / 2)) / np.cos(grid_latlon_radians[:, 0]) 
# Compute scale factor for each grid cell
k0 = pstere_lat2k0(srs.GetProjParm('latitude_of_origin'))
scale_factors = np.sqrt((m * k0))
# reshape scale_factors to the same dimensions as the grid
scale_factors = scale_factors.reshape(xi.shape)

# Create lat/lon variables of projected grid
grid_lon_2d = grid_lon.reshape(xi.shape)
grid_lat_2d = grid_lat.reshape(yi.shape)

## Begin the interpolation/projection

# Are we just testing the first year or not?
if not args.test:
	process_times = times
else:
	date_end_here = (dt.datetime.strptime(args.date_start, '%Y-%m-%d') + dt.timedelta(days=364)).strftime('%Y-%m-%d')
	process_times = pd.date_range(args.date_start, date_end_here, freq='1MS')

# Set up store in which to hold projected (gridded) data
store = np.zeros((len(process_times), yi.shape[0], yi.shape[1]))
# Integer time counter for indexing to np store array
nt = 0

# Do the projection to grid
for t in process_times:

	print(t.strftime('%Y-%m'))

	# Extract runoff mmWE value at each coordinate
	r = runoff.runoff.sel(time=t).values
	# And flatten to a vector which matches the xy (grid coords) array order
	r = r.flatten()

	# Interpolate point data onto grid
	zi = interpolate.griddata(xy, r, (xi, yi), method='linear')

	# Apply scale_factor
	zi_sf = zi / scale_factors

	# Save to store array
	store[nt, :, :] = zi_sf

	# Increment year counter
	nt += 1


## Create xarray data array
# Coordinates of the data, for xarray representation
coords = {'TIME': process_times, 'Y': yi[:,0], 'X': xi[0,:]}
# Projected runoff data
da = xr.DataArray(store, coords=coords, 
	dims=['TIME', 'Y', 'X'], encoding={'dtype':np.dtype('Float32')})

# Set attributes from RACMO netCDF
da.name = runoff.runoff.long_name
da.attrs['long_name'] = runoff.runoff.long_name
da.attrs['units'] = runoff.runoff.units
da.attrs['standard_name'] = runoff.runoff.standard_name
da.attrs['grid_mapping'] = 'polar_stereographic'


## Create associated lat/lon coordinates DataArrays
coords_geo = {'Y': yi[:,0], 'X': xi[0,:]}

lon_da = xr.DataArray(grid_lon_2d, coords=coords_geo, dims=['Y', 'X'], 
	encoding={'_FillValue': -9999.})
lon_da.attrs['grid_mapping'] = 'polar_stereographic'
lon_da.attrs['units'] = 'degrees'
lon_da.attrs['standard_name'] = 'longitude'

lat_da = xr.DataArray(grid_lat_2d, coords=coords_geo, dims=['Y', 'X'], 
	encoding={'_FillValue': -9999.})
lat_da.attrs['grid_mapping'] = 'polar_stereographic'
lat_da.attrs['units'] = 'degrees'
lat_da.attrs['standard_name'] = 'latitude'


## Define projection
crs = xr.DataArray(0, encoding={'dtype':np.dtype('int8')})
crs.attrs['grid_mapping_name'] = 'polar_stereographic'
crs.attrs['scale_factor_at_central_origin'] = srs.GetProjParm('scale_factor')
crs.attrs['standard_parallel'] = srs.GetProjParm('latitude_of_origin')
crs.attrs['straight_vertical_longitude_from_pole'] = srs.GetProjParm('central_meridian')
crs.attrs['false_easting'] = srs.GetProjParm('false_easting')
crs.attrs['false_northing'] = srs.GetProjParm('false_northing')
# lat_0 does not have a WKT representation for polar stereographic
# http://cfconventions.org/wkt-proj-4.html
# look it up from Proj4 string directly instead
p4 = srs.ExportToProj4()
lat_0_pos = p4.find('+lat_0=')
lat_0 = float(p4[lat_0_pos+7:lat_0_pos+7+3].strip())
crs.attrs['latitude_of_projection_origin'] = lat_0


## Create the Dataset
ds = xr.Dataset({'runoff':da, 'lon': lon_da, 'lat': lat_da, 'polar_stereographic':crs})

# Main metadata
ds.attrs['Conventions'] = 'CF-1.4'
ds.attrs['history'] = 'This NetCDF generated using bitbucket atedstone/fwf/reproject_racmo.py on %s provided by JLB/BN/MvdB\n %s' %(args.fn_RACMO, runoff.history)
ds.attrs['institution'] = 'University of Bristol (Andrew Tedstone), IMAU (Brice Noel)'
ds.attrs['title'] = 'Monthly runoff in the RACMO 2.3 domain on a projected grid'
ds.attrs['source'] = 'RACMO 2.3'
ds.attrs['proj4'] = grid_proj.srs

# Additional geo-referencing
ds.attrs['nx'] = float(xi.shape[1])
ds.attrs['ny'] = float(yi.shape[0])
ds.attrs['xmin'] = float(np.round(np.min(xi), 0))
ds.attrs['ymax'] = float(np.round(np.max(yi), 0))
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


## Save the Dataset as a NetCDF file
if args.test:
	fn_save = args.fn_RACMO.split('/')[-1][:-3] + '_%s_TEST.nc' % args.grid
else:
	fn_save = args.fn_RACMO.split('/')[-1][:-3] + '_%s.nc' % args.grid
print(PROCESS_DIR + fn_save)
ds.to_netcdf(PROCESS_DIR + fn_save,	format='NETCDF4')



##############################################################################
## Do some simple checks if projecting the bamber grid (as per used in 2012 paper)

if args.grid == 'bamber':
	# To create this file, use Williams mc_land_mask___bamber_proj.tif:
	# gdalwarp -tr 5000 5000 -te -800000 -3400000 700000 -600000 mc_land_mask___bamber_proj.tif mc_land_mask___bamber_proj_5km.tif
	mask = georaster.SingleBandRaster('/scratch/bedmachine/mc_land_mask___bamber_proj_5km.tif')
	ice_area = np.where(mask.r == 2, True, False)

	# Bamber's value for gridded product stated as 242.888 km^3 in Interpolate_racmo.BAK
	print('1958 runoff flux:')
	print(((ds.runoff.sel(TIME=slice('1958-01-01','1958-12-01')) * np.flipud(ice_area)) * (5*5) / 1.0e6).sum())



##############################################################################
## Project masks(s) etc to geotiff
# Use all the same grid logic as defined above
print('Projecting mask(s)...')

# Set geoTransform for writing to geotiff
trans = (xi[0,0], (xi[0,1]-xi[0,0]), 0, yi[-1,0], 0, (yi[-2,0]-yi[-1,0]))

# Set up coords for saving to netcdf
masks_coords = {'Y': yi[:,0], 'X': xi[0,:]}

# Open RACMO masks (they are on same grid)
masks_ds = xr.open_dataset(args.fn_RACMO_masks)


## Mask of all ice areas; we'll need to split into numerical basins later
mask_ice_pts = masks_ds.icecon.values
# And flatten to a vector which matches the xy (grid coords) array order
mask_ice_pts = mask_ice_pts.flatten()

# Interpolate point data onto grid
mask_ice_gridded = interpolate.griddata(xy, mask_ice_pts, (xi, yi), method='nearest')

# Write to geotiff
georaster.simple_write_geotiff(
	fn_mask_ice,
	np.flipud(mask_ice_gridded), 
	trans,
	proj4=grid_proj.srs,
	dtype=gdal.GDT_Byte
	)
xr_mask_ice = xr.DataArray(mask_ice_gridded, dims=('Y', 'X'), coords=masks_coords)


## Ice-only Topography
topo_ice_pts = masks_ds.topography.where(masks_ds.icecon == 1).values
# And flatten to a vector which matches the xy (grid coords) array order
topo_ice_pts = topo_ice_pts.flatten()

# Interpolate point data onto grid
topo_ice_gridded = interpolate.griddata(xy, topo_ice_pts, (xi, yi), method='linear')

# Write to geotiff
georaster.simple_write_geotiff(
	fn_dem_ice,
	np.flipud(topo_ice_gridded), 
	trans,
	proj4=grid_proj.srs,
	dtype=gdal.GDT_Float64
	)
xr_topo_ice = xr.DataArray(topo_ice_gridded, dims=['Y','X'], coords=masks_coords, encoding={'dtype':np.dtype('Float32')})


## Land-and-ice mask
mask_landandice_pts = (masks_ds.LSM_noGrIS + masks_ds.Gr_land).values
mask_landandice_pts = mask_landandice_pts.flatten()
mask_landandice_gridded = interpolate.griddata(xy, mask_landandice_pts, (xi, yi), method='nearest')
mask_landandice_gridded = ndimage.binary_fill_holes(mask_landandice_gridded)

# Write to geotiff
georaster.simple_write_geotiff(
	fn_mask_landandice,
	np.flipud(mask_landandice_gridded),
	trans,
	proj4=grid_proj.srs,
	dtype=gdal.GDT_UInt16
	)
xr_mask_landandice = xr.DataArray(mask_landandice_gridded, dims=['Y','X'], coords=masks_coords)


## Land-only mask
mask_land_pts = ((masks_ds.LSM_noGrIS + masks_ds.Gr_land) - (masks_ds.icecon == 1)).values
mask_land_pts = mask_land_pts.flatten()
mask_land_gridded = interpolate.griddata(xy, mask_land_pts, (xi, yi), method='nearest')

# Write to geotiff
georaster.simple_write_geotiff(
	fn_mask_land,
	np.flipud(mask_land_gridded),
	trans,
	proj4=grid_proj.srs,
	dtype=gdal.GDT_UInt16
	)
xr_mask_land = xr.DataArray(mask_land_gridded, dims=['Y','X'], coords=masks_coords)


## Greenland land+ice
greenland_pts = masks_ds.Gr_land.values
greenland_pts = greenland_pts.flatten()
greenland_gridded = interpolate.griddata(xy, greenland_pts, (xi, yi), method='nearest')
greenland_gridded = np.where(greenland_gridded > 0, 1, 0)
# Write to geotiff
georaster.simple_write_geotiff(
	fn_mask_greenland,
	np.flipud(greenland_gridded),
	trans,
	proj4=grid_proj.srs,
	dtype=gdal.GDT_UInt16
	)
xr_mask_greenland = xr.DataArray(greenland_gridded, dims=['Y','X'], coords=masks_coords)



## The others for netcdf only
icemask_noGrIS = masks_ds.icemask_no_GrIS.values
icemask_noGrIS_gridded = interpolate.griddata(xy, icemask_noGrIS.flatten(), (xi, yi), method='nearest')
xr_mask_icemask_noGrIS = xr.DataArray(icemask_noGrIS_gridded, dims=['Y','X'], coords=masks_coords)

LSM_noGrIS = masks_ds.LSM_noGrIS.values
LSM_noGrIS_gridded = interpolate.griddata(xy, LSM_noGrIS.flatten(), (xi, yi), method='nearest')
xr_mask_LSM_noGrIS = xr.DataArray(LSM_noGrIS_gridded, dims=['Y','X'], coords=masks_coords)

Gr_land = masks_ds.Gr_land.values
Gr_land_gridded = interpolate.griddata(xy, Gr_land.flatten(), (xi, yi), method='nearest')
xr_mask_Gr_land = xr.DataArray(Gr_land_gridded, dims=['Y','X'], coords=masks_coords)

GrIS_mask = masks_ds.GrIS_mask.values
GrIS_mask_gridded = interpolate.griddata(xy, GrIS_mask.flatten(), (xi, yi), method='nearest')
xr_mask_GrIS_mask = xr.DataArray(GrIS_mask_gridded, dims=['Y','X'], coords=masks_coords)

GrIS_caps_mask = masks_ds.GrIS_caps_mask.values
GrIS_caps_mask_gridded = interpolate.griddata(xy, GrIS_caps_mask.flatten(), (xi, yi), method='nearest')
xr_mask_GrIS_caps_mask = xr.DataArray(GrIS_caps_mask_gridded, dims=['Y','X'], coords=masks_coords)

## Also put them all into a netCDF
new_ds_masks = xr.Dataset({
	'LSM_noGrIS':xr_mask_LSM_noGrIS,
	'icemask_no_GrIS':xr_mask_icemask_noGrIS,
	'Gr_land':xr_mask_Gr_land,
	'GrIS_mask':xr_mask_GrIS_mask,
	'GrIS_caps_mask':xr_mask_GrIS_caps_mask,
	'Greenland_all':xr_mask_greenland,
	'all_land_no_ice':xr_mask_land,
	'landandice':xr_mask_landandice,
	'all_ice':xr_mask_ice,
	'ice_topo':xr_topo_ice
	})

fn_save = args.fn_RACMO_masks.split('/')[-1][:-3] + '_%s.nc' % args.grid
new_ds_masks.to_netcdf(PROCESS_DIR + fn_save)


print('Finished.')


"""
Code to plot all masks:

for i in masks.variables:
	print(i)
	if i in ('X','Y'): continue
	figure(),masks[i].plot(),title(i)

"""

#############################################################################

## An example runoff file for testing runoff routing
# runoff_195807 = ds.runoff.sel(TIME='1958-07-01').where(mask_gridded == 1)
# runoff_195807 = runoff_195807 * 5*5 / 1.0e6
# georaster.simple_write_geotiff(
# 	'/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/racmo_pstere_5km_runoff_1958-07.tif',
# 	np.flipud(runoff_195807),
# 	trans,
# 	proj4=grid_proj.srs,
# 	dtype=gdal.GDT_Float64
# 	)
