"""
Project RACMO 2.3 data onto a regular grid using interpolation.
Outputs an equivalent NetCDF of the time series.

Remember to set environment variable PROCESS_DIR

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
from osgeo import osr
from scipy import ndimage
import os
import argparse

# obtain from https://github.com/atedstone/georaster
import georaster

parser = argparse.ArgumentParser(description='Project geographic RACMO data for Arctic domain')

# Positional arguments
parser.add_argument('fn_RACMO', type=str, help='str, absolute path to RACMO netCDF file')
parser.add_argument('fn_RACMO_masks', type=str, help='str, absolute path to RACMO masks netCDF file')
parser.add_argument('date_start', type=str, help='str, start date of RACMO data in format yyyy-mm-dd')
parser.add_argument('date_end', type=str, help='str, end date of RACMO data in format yyyy-mm-dd')

parser.add_argument('-grid', type=str, dest='grid', default='pstere', help='str, grid to use, pstere or bamber')
parser.add_argument('-test', type=bool, dest='test', action='store_true', help='Test script (run only for first year of data)')

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
if args.grid == 'pstere':
	grid_proj = pyproj.Proj('+init=EPSG:3413')
	fn_mask_ice = PROCESS_DIR + 'mask_ice_racmo_EPSG3413_5km.tif'
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


## Create grid at 5 km resolution
if grid == 'pstere':
	yi, xi = np.mgrid[np.min(y):np.max(y):5000, np.min(x):np.max(x):5000]
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


## Begin the interpolation/projection
# Set up store in which to hold projected (gridded) data
store = np.zeros((len(times), yi.shape[0], yi.shape[1]))
# Integer time counter for indexing to np store array
nt = 0

# Are we just testing the first year or not?
if not args.test:
	process_times = times
else:
	date_end_here = (dt.strptime(args.date_start, '%Y-%m-%d') + dt.timedelta(days=365).strftime('%Y-%m-%d'))
	process_times = pd.date_range(args.date_start, date_end_here, freq='1MS')

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
coords = {'TIME': times, 'Y': yi[:,0], 'X': xi[0,:]}
da = xr.DataArray(store, coords=coords, dims=['TIME', 'Y', 'X'])
# Set attributes from RACMO netCDF
da.name = runoff.runoff.long_name
da.attrs['long_name'] = runoff.runoff.long_name
da.attrs['units'] = runoff.runoff.units
da.attrs['standard_name'] = runoff.runoff.standard_name


## Create xarray dataset (for saving to netCDF)
ds = xr.Dataset({'runoff':da})
# Include history from RACMO netCDF
ds.attrs['history'] = runoff.history + ' // This NetCDF generated using bitbucket atedstone/fwf/reproject_racmo.py on %s provided by JLB/BN/MvdB' % args.fn_RACMO
ds.attrs['author'] = 'Andrew Tedstone a.j.tedstone@bristol.ac.uk'
# ds.attrs['GeoTransform'] = (xi[0,0], (xi[0,1]-xi[0,0]), 0, yi[-1,0], 0, (yi[-2,0]-yi[-1,0]))
# ds.attrs['Conventions'] = "CF-1.6"

# ds.X.attrs['units'] = 'meters'
# ds.X.attrs['long_name'] = 'Eastings'
# ds.Y.attrs['units'] = 'meters'
# ds.Y.attrs['long_name'] = 'Northings'

if args.test:
	fn_save = args.fn_RACMO[:-3] + '_%s_TEST.nc' % args.grid
else:
	fn_save = args.fn_RACMO[:-3] + '_%s.nc' % args.grid
ds.to_netcdf(PROCESS_DIR + fn_save,	format='NETCDF4')



##############################################################################
## Do some simple checks if projecting the bamber grid (as per used in 2012 paper)

if args.grid == 'bamber':
	mask = georaster.SingleBandRaster('/media/sf_Data/mc_land_mask___bamber_proj_5km.tif')
	ice_area = np.where(mask.r == 2, True, False)

	# Bamber's value for gridded product stated as 242.888 km^3 in Interpolate_racmo.BAK
	print('1958 runoff flux:')
	print(((ds.runoff.sel(TIME=slice('1959-01-01','1959-09-01')) * ice_area) * (5*5) / 1.0e6).sum())



##############################################################################
## Project masks(s) etc to geotiff
# Use all the same grid logic as defined above
print('Projecting mask(s)...')

# Set geoTransform for writing to geotiff
trans = (xi[0,0], (xi[0,1]-xi[0,0]), 0, yi[-1,0], 0, (yi[-2,0]-yi[-1,0]))

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
	dtype=georaster.gdal.GDT_Byte
	)


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
	dtype=georaster.gdal.GDT_Float64
	)


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
	dtype=georaster.gdal.GDT_UInt16
	)

print('Finished.')

#############################################################################

## An example runoff file for testing runoff routing
# runoff_195807 = ds.runoff.sel(TIME='1958-07-01').where(mask_gridded == 1)
# runoff_195807 = runoff_195807 * 5*5 / 1.0e6
# georaster.simple_write_geotiff(
# 	'/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/racmo_pstere_5km_runoff_1958-07.tif',
# 	np.flipud(runoff_195807),
# 	trans,
# 	proj4=grid_proj.srs,
# 	dtype=georaster.gdal.GDT_Float64
# 	)
