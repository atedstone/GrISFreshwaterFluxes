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
from osgeo import gdal, osr
import math
from scipy.ndimage import morphology as sc_morph
import subprocess

import georaster

# If true, only process 1958
debug = True

PROCESS_DIR = os.environ['PROCESS_DIR']


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

# ----------------------------------------------------------------------------

## Only use topography of ice surfaces - route to ice sheet + cap margins.
dem_uri = PROCESS_DIR + 'mask_Geopotential_mosaic.tif'
ice_mask_uri = PROCESS_DIR + 'mask_icemask_mosaic.tif'
ice_mask = georaster.SingleBandRaster(ice_mask_uri)
# Convert PROMICE indicators to binary
ice_mask.r = np.where(ice_mask.r >= 1, 1, 0)
# Fill the small holes which arise from the downsampling
ice_mask_filled = sc_morph.binary_fill_holes(ice_mask.r)
#filled_areas = filled_mask - ice_mask.r
#ice_mask.r = filled_mask
ice_mask.save_geotiff(PROCESS_DIR + 'mask_icemask_mosaic_filled_binary.tif')
dem_sm_uri = PROCESS_DIR + 'topography_ice_mosaic_smoothed.tif'
## Set nans to zero - I think this is needed for pit filling?
dem = georaster.SingleBandRaster(dem_uri)
dem.r = np.where(dem.r == 9999, 0, dem.r)
# Add 200 m elevation to areas of ice mask which got filled in
#dem.r = np.where(filled_areas, dem.r + 200, dem.r)
dem_ice = georaster.SingleBandRaster.from_array(np.where(ice_mask.r == 1, dem.r, 0), ice_mask.trans,
	ice_mask.proj.srs, gdal.GDT_Float32)
#dem.r = np.where(np.isnan(dem.r), 0, dem.r)
dem_ice.save_geotiff(dem_sm_uri)
dem = None


## Step 0: fill pits
dempits_uri = PROCESS_DIR + 'ROUTING_racmo_EPSG3413_5km_dem_filledpits.tif'
routing.fill_pits(dem_sm_uri, dempits_uri)

# Generate revised ice mask as pit-filling changes margins
# dempits = georaster.SingleBandRaster(dempits_uri)
# ice_mask_pits = georaster.SingleBandRaster.from_array(
# 	np.where(dempits.r > 0, 1, 0),
# 	ice_mask.trans, ice_mask.proj.srs, gdal.GDT_Float32
# 	)
# # Fill holes again...seems to be necessary to deal with fjord artifacts on east coast
# ice_mask_pits.r = sc_morph.binary_fill_holes(ice_mask_pits.r)
# ice_mask_pits.save_geotiff(PROCESS_DIR + 'mask_icemask_mosaic_filled_binary_pits.tif')


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
landice_mask = georaster.SingleBandRaster(PROCESS_DIR + 'mask_LandSeaMask.tif')
trans = landice_mask.trans
proj4 = landice_mask.srs.ExportToProj4()
# Also output filled Greenland-only mask (faciliates easy comparison to unrouted data)
Gr_land = georaster.SingleBandRaster(PROCESS_DIR + 'mask_LSMGr.tif')
Gr_land_filled_r = sc_morph.binary_fill_holes(Gr_land.r)
Gr_land_filled = georaster.SingleBandRaster.from_array(Gr_land_filled_r, trans,
	landice_mask.proj.srs, gdal.GDT_Byte)
Gr_land_filled.save_geotiff(PROCESS_DIR + 'mask_LSMGr_filled.tif',
	dtype=gdal.GDT_Byte)

#landice_mask.r = np.where(((landice_mask.r ==1) | (Gr_land_filled.r == 1)), 1, 0)
# Fill holes to make sure we reach the coast
landice_mask_filled = sc_morph.binary_fill_holes(landice_mask.r)
landice_mask.r = landice_mask_filled
landice_mask.save_geotiff(PROCESS_DIR + 'mask_LandSeaMask_filled.tif',
	dtype=gdal.GDT_Byte)
landice_mask = None
# Merge in the wider extent of Gr_land_filled
subprocess.call('gdalwarp -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/gris_bounds.shp %s/mask_LSMGr_filled.tif %s/mask_LandSeaMask_filled.tif' %(PROCESS_DIR, PROCESS_DIR), shell=True)
# Load the filled one properly
landice_mask_filled = georaster.SingleBandRaster(PROCESS_DIR + 'mask_LandSeaMask_filled.tif')

# landice_mask = None
# landice_mask_filled = None
# # LSMGr and LandSeaMask have different Greenland land extents
# # Overwrite LandSeaMask with LSMGr so that we can isolate GrIS in post-processing
# subprocess.call('gdalwarp -r max %s/mask_LSMGr_filled.tif %s/mask_LandSeaMask_filled.tif' %(PROCESS_DIR, PROCESS_DIR), shell=True)
# landice_mask_filled = georaster.SingleBandRaster(PROCESS_DIR + 'mask_LandSeaMask_filled.tif')
# Get distances
out, distances = morphology.medial_axis(landice_mask_filled.r, return_distance=True)
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

## Load land, no ice mask to route tundra runoff (euclidean distance)
# land_mask = georaster.SingleBandRaster.from_array(landice_mask_filled.r-ice_mask_filled, trans,
# 	proj4, gdal.GDT_Byte)
# land_mask.save_geotiff(PROCESS_DIR + 'ROUTING_land_only_mask.tif')
#icemask_11km = georaster.SingleBandRaster(PROCESS_DIR + 'mask_icemask.tif')
landmask_11km = georaster.SingleBandRaster(PROCESS_DIR + 'mask_LandSeaMask.tif')
# Use the 'detailed' (derived from 1km) ice mask to identify just-land pixels.
# Note that this isn' the filled version, because the ice FWF routing also isn't over filled areas.
# However, am using 11 km land mask as nothing better is available from IMAU datasets provided.
land_only_11km = landmask_11km.r - ice_mask.r
# x = np.arange(0, land_mask.nx)
# y = np.arange(0, land_mask.ny)
x = np.arange(0, icemask_11km.nx)
y = np.arange(0, icemask_11km.ny)
xi, yi = np.meshgrid(x, y)
xp = xi[land_only_11km == 1].flatten()
yp = yi[land_only_11km == 1].flatten()
land_points = np.zeros((len(xp), 2))
land_points[:, 0] = yp
land_points[:, 1] = xp


### Calculate scaling factors for polar stereo projection
# We also use these grids later as nc variables lon and lat
grid_lon, grid_lat = landice_mask_filled.coordinates(latlon=True)
# Convert to paired columns
grid_latlon = np.vstack((grid_lat.flatten(), grid_lon.flatten())).T
# Convert to radians
grid_latlon_radians = grid_latlon / 57.29578
# Calculate latitudinal scaling
m = 2.0 * np.tan(45.0 / 57.29578 - (grid_latlon_radians[:, 0] / 2)) / np.cos(grid_latlon_radians[:, 0]) 
# Compute scale factor for each grid cell
k0 = pstere_lat2k0(70.)
scale_factors = np.power(m * k0, 2)
# reshape scale_factors to the same dimensions as the grid
scale_factors = scale_factors.reshape(grid_lon.shape)


## Step 3: route monthly runoff
times = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')

if debug == True:
	store_ice = np.zeros((12*1, ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	store_tundra = np.zeros((12*1, ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	dates = pd.date_range('1958-01-01', '1958-12-01', freq='1MS')
else:
	dates = times
	store_ice = np.zeros((len(dates), ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))
	store_tundra = np.zeros((len(dates), ice_mask.ds.RasterYSize, ice_mask.ds.RasterXSize))

n = 0
ice_r_pre = []
ice_r_post = []
for date in dates:
	print(date)
	# Import runoff from reprojected GeoTIFF
	fname = PROCESS_DIR + 'runoff_b%s.tif' % (n+1)
	r_month = georaster.SingleBandRaster(fname).r 
	# Convert mmWE to flux per grid cell
	r_month = np.where(r_month > 0, r_month * (5*5) / 1.0e6, 0)  / scale_factors 
	# Get ice component of runoff (USING FILLED MASK! - IMPORTANT)
	r_month_ice = np.where(ice_mask_filled == 1, r_month, np.nan)
	# Calculate runoff value pre-routing for debugging
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
	# store_ice[n,:,:] = np.flipud(coast_grid_ice)
	store_ice[n,:,:] = coast_grid_ice

	# MUST close handle to flux file otherwise pygeoprocessing can't overwrite on next iteration
	flux = None


	## Route tundra fluxes too ...
	## NEED TO OPEN UNMODIFIED 11KM FILES HERE
	# Do by euclidean distance
	r_base11km = georaster.SingleBandRaster(PROCESS_DIR + 'arctic_domain_bckup/' + 'runoff_b%s.tif' %(n+1))
	r_base11km.r = (r_base11km.r * (5*5) / 1.0e6) / scale_factors 
	coast_grid_tundra = np.zeros((dist_land.ny, dist_land.nx))
	for ry, rx in land_points:
		distance, index = tree.query((ry, rx), k=1)
		cpy, cpx = coast_points[index, :]
		coast_grid_tundra[int(cpy), int(cpx)] += r_base11km.r[int(ry), int(rx)]

	r_base11km = None
	# Save coastal fluxes to store
	store_tundra[n,:,:] = coast_grid_tundra

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

grid_x,grid_y = dem_ice.coordinates()
coords = {'TIME':dates, 'Y':grid_y[:,0], 'X':grid_x[0,:]}

da_ice = xr.DataArray(np.round(store_ice, 2), 
	coords=coords, 
	dims=['TIME', 'Y', 'X'], 
	encoding={'dtype':'int16', 'scale_factor':0.01, 'zlib':True, '_FillValue':-9999})
da_ice.name = 'Ice sheet runoff'
da_ice.attrs['long_name'] = 'Ice sheet runoff'
da_ice.attrs['units'] = 'km3'
da_ice.attrs['grid_mapping'] = 'polar_stereographic'

da_tundra = xr.DataArray(np.round(store_tundra, 2), 
	coords=coords, 
	dims=['TIME', 'Y', 'X'], 
	encoding={'dtype':'int16', 'scale_factor':0.01, 'zlib':True, '_FillValue':-9999})
da_tundra.name = 'Tundra runoff'
da_tundra.attrs['long_name'] = 'Tundra runoff'
da_tundra.attrs['units'] = 'km3'
da_tundra.attrs['grid_mapping'] = 'polar_stereographic'


## Define projection
srs = osr.SpatialReference()
srs.ImportFromProj4('+init=epsg:3413')
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

## Create associated lat/lon coordinates DataArrays
coords_geo = {'Y': coords['Y'], 'X': coords['X']}

lon_da = xr.DataArray(grid_lon, coords=coords_geo, dims=['Y', 'X'], 
	encoding={'_FillValue': -9999.})
lon_da.attrs['grid_mapping'] = 'polar_stereographic'
lon_da.attrs['units'] = 'degrees'
lon_da.attrs['standard_name'] = 'longitude'

lat_da = xr.DataArray(grid_lat, coords=coords_geo, dims=['Y', 'X'], 
	encoding={'_FillValue': -9999.})
lat_da.attrs['grid_mapping'] = 'polar_stereographic'
lat_da.attrs['units'] = 'degrees'
lat_da.attrs['standard_name'] = 'latitude'


ds = xr.Dataset({'runoff_ice':da_ice, 
	'runoff_tundra':da_tundra,
	'lon':lon_da,
	'lat':lat_da,
	'polar_stereographic':crs})

# Main metadata
ds.attrs['Conventions'] = 'CF-1.4'
ds.attrs['history'] = 'This NetCDF generated using bitbucket atedstone/fwf/arctic_runoff_routing.py using data output by fwf/project_racmo.R on /scratch/L0data/RACMO/RACMO2.3_GRN11_runoff_monthly_1958-2015.nc'
ds.attrs['institution'] = 'University of Bristol (Andrew Tedstone)'
ds.attrs['title'] = 'Monthly ice sheet and tundra runoff routed to coastal pixels'

# Additional geo-referencing
ds.attrs['nx'] = float(dempits.nx)
ds.attrs['ny'] = float(dempits.ny)
ds.attrs['xmin'] = float(np.round(np.min(grid_x), 0))
ds.attrs['ymax'] = float(np.round(np.max(grid_y), 0))
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

ds.to_netcdf('/home/at15963/Dropbox/work/papers/bamber_fwf/outputs_Nov2017/FWF17_runoff_RACMO2.3p2.nc', format='NetCDF4')

pre = pd.Series(ice_r_pre, index=dates)
pre.to_csv(PROCESS_DIR + 'runoff_monthly_totals_pre_routing.csv')
post = pd.Series(ice_r_post, index=dates)
post.to_csv(PROCESS_DIR + 'runoff_monthly_totals_post_routing.csv')
