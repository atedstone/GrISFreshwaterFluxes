"""

Allocate coastal outflux pixels to ocean basins

Uses IHO World Seas definitions.

@author Andrew Tedstone a.j.tedstone@bristol.ac.uk
@created 2017-February-07

"""

import numpy as np
import geopandas as gpd
import georaster
from shapely.geometry import Point, Polygon
import os

PROCESS_DIR = os.environ['PROCESS_DIR']

dist_land = georaster.SingleBandRaster(PROCESS_DIR + 'ROUTING_distances_landandice.tif')
coast = np.where(dist_land.r == 1, 1, 0)

iho = gpd.read_file('/home/at15963/Dropbox/work/data/World_Seas_IHO_v2/World_Seas_IHO_v2.shp')
# Retain only the basins of interest - makes the distance calculations later faster
in_arctic = iho.intersects(Polygon(((-121, 52), (43, 52), (43, 90), (-121, 90))))
iho_arctic = iho[in_arctic]
iho = None

# Get grids of lat/lon coordinates 
lons, lats = dist_land.coordinates(latlon=True)
lons_coast = np.where(coast == 1, lons, np.nan)
lats_coast = np.where(coast == 1, lats, np.nan)
lons_coast_pts = lons_coast[~np.isnan(lons_coast)].flatten()
lats_coast_pts = lats_coast[~np.isnan(lats_coast)].flatten()
ll_coast_points = np.zeros((len(lons_coast_pts), 2))
ll_coast_points[:, 0] = lons_coast_pts
ll_coast_points[:, 1] = lats_coast_pts

# Calculate equivalent pixel coordinates so we can create a raster of outflux basins
x = np.arange(0, dist_land.nx)
y = np.arange(0, dist_land.ny)
xi, yi = np.meshgrid(x, y)
xp = xi[np.where(coast == 1, True, False)].flatten()
yp = yi[np.where(coast == 1, True, False)].flatten()
coast_points = np.zeros((len(xp), 2))
coast_points[:, 0] = yp
coast_points[:, 1] = xp

basins = []
basins_raster = np.zeros((dist_land.ny, dist_land.nx))
for (y, x), (lon, lat) in zip(coast_points, ll_coast_points):

	p = Point(lon, lat)
	distances = iho_arctic.distance(p)
	idx = distances.idxmin()
	basins.append(idx)
	basins_raster[y, x] = idx

basin_desc = ''
for b in iho_arctic.iterrows():
	basin_desc += ', %s:%s' % (b[0], b[1]['NAME'])

history = 'Basins are as defined by "Limits of Oceans & Seas, Special Publication No. 23" published by the IHO in 1953. The dataset was composed by the Flanders Marine Data and Information Centre and downloaded from www.marineregions.org on 6 February 2017. Coastal runoff pixels were allocated an ocean basin to run off into using fwf/allocate_coast_to_basin.py.'
author = 'Andrew Tedstone, University of Bristol'
metadata = { 'basins':basin_desc, 'history':history, 'author':author }

## Save basins raster
georaster.simple_write_geotiff(
	'/scratch/process/outflow_basins.tif',
	basins_raster, 
	dist_land.trans,
	proj4=dist_land.srs.ExportToProj4(),
	dtype=georaster.gdal.GDT_UInt16,
	metadata=metadata
	)


