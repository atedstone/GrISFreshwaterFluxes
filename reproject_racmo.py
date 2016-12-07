"""
Project RACMO 2.3 data onto a regular grid using interpolation.
Outputs an equivalent NetCDF of the time series.
Filenames are hard-coded into the script.

@author Andrew Tedstone a.j.tedstone@bristol.ac.uk
@created 2016-December-05

"""

import xarray as xr
import pyproj
from scipy import interpolate
import pandas as pd
import numpy as np
from osgeo import osr

# obtain from https://github.com/atedstone/georaster
import georaster

# Set grid to 'pstere' or 'bamber'
# 'bamber' grid is GrIS only, 'pstere' will result in grid across full RACMO extent
# Use 'bamber' grid only for testing to ensure that interpolation result matches 2012 FWF paper results
# Characteristics of each grid are hard-coded into logic below
grid = 'pstere'

# Run interpolation for complete time series (TRUE), or just first year (1958, FALSE)?
# Use for testing 
complete_ts = False


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
if grid == 'pstere':
	grid_proj = pyproj.Proj('+init=EPSG:3413')
	mask_out_file = '/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/mask_racmo_EPSG3413_5km.tif'
elif grid == 'bamber':
	grid_proj = pyproj.Proj('+proj=sterea +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
	mask_out_file = '/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/mask_racmo_bamberpstere_5km.tif'

# Create SRS representation
srs = osr.SpatialReference()
srs.ImportFromProj4(grid_proj.srs)


## Open RACMO runoff dataset
runoff = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_runoff_monthly_1958-2015.nc',
		decode_times=False, chunks={'time':15})
# Apply time dimension in format understandable by xarray
times = pd.date_range('1958-01-01', '2015-12-01', freq='1MS')
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
print('k0 calculated as %.5f' % k0)
scale_factors = np.sqrt((m * k0))
# reshape scale_factors to the same dimensions as the grid
scale_factors = scale_factors.reshape(xi.shape)


## Begin the interpolation/projection
# Set up store in which to hold projected (gridded) data
store = np.zeros((len(times), yi.shape[0], yi.shape[1]))
# Integer time counter for indexing to np store array
nt = 0

# Are we just testing the first year or not?
if complete_ts:
	process_times = times
else:
	process_times = pd.date_range('1958-01-01', '1958-12-01', freq='1MS')

for t in process_times:

	print(t.strftime('%Y-%m'))

	# Extract runoff mmWE value at each coordinate
	r = runoff.runoff.sel(time=t).values
	# And flatten to a vector which matches the xy (grid coords) array order
	r = r.flatten()

	# Interpolate point data onto grid
	zi = interpolate.griddata(xy, r, (xi, yi), method='cubic')

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
ds.attrs['history'] = runoff.history + ' // This NetCDF generated using bitbucket atedstone/fwf/reproject_racmo.py on RACMO2.3_GRN11_runoff_monthly_1958-2015.nc provided by JLB/BN/MvdB'
ds.attrs['author'] = 'Andrew Tedstone a.j.tedstone@bristol.ac.uk'
ds.attrs['gridding method'] = 'cubic spline scipy.interpolate.griddata'
ds.to_netcdf('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_EPSG3413_runoff_monthly_1958.nc')



##############################################################################
## Do some simple checks if projecting the bamber grid (as per used in 2012 paper)

if grid == 'bamber':
	mask = georaster.SingleBandRaster('/media/sf_Data/mc_land_mask___bamber_proj_5km.tif')
	ice_area = np.where(mask.r == 2, True, False)

	# Bamber's value for gridded product stated as 242.888 km^3 in Interpolate_racmo.BAK
	print('1958 runoff flux:')
	print(((ds.runoff.sel(TIME=slice('1959-01-01','1959-09-01')) * ice_area) * (5*5) / 1.0e6).sum())



##############################################################################
## Project masks(s)
# Use all the same grid logic as defined above
print('Projecting mask(s)...')

## Open RACMO masks (they are on same grid)
masks_ds = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_masks.nc')

# Select mask of all ice areas; we'll need to split into numerical basins later
mask_racmo = masks_ds.icecon.values
# And flatten to a vector which matches the xy (grid coords) array order
mask_racmo = mask_racmo.flatten()

# Interpolate point data onto grid
zi = interpolate.griddata(xy, mask_racmo, (xi, yi), method='nearest')

# Write to geotiff
trans = (xi[0,0], (xi[0,1]-xi[0,0]), 0, yi[0,0], 0, (yi[1,0]-yi[0,0]))
georaster.simple_write_geotiff(
	mask_out_file,
	zi, 
	trans,
	proj4=grid_proj.srs,
	dtype=georaster.gdal.GDT_Byte
	)

