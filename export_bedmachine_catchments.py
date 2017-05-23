"""
Convert bedmachine named points (provided by Chris Williams, May 2017) to
shapefile of polygons, which attached name, to help with identifying outflow
basins
"""

import pandas as pd
import geopandas as gpd 
from shapely.geometry import Point, MultiPoint

data = pd.read_csv('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/bedmachine_all_points_named.csv',
	names=('x', 'y', 'n1', 'n2', 'name'), usecols=('x', 'y', 'name'))
data.name = data.name.str.strip()

# These points are in Morlighem projection
proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

geometry = [Point(xy) for xy in zip(data.x, data.y)]
data_geo = gpd.GeoDataFrame(data.name, geometry=geometry, crs=proj4)

data_latlon = data_geo.to_crs({'init':'epsg:4326'})

# Do a groupby, which is a plain Pandas operation...
multipoint = data_latlon.groupby(data_latlon.name).apply(lambda x: MultiPoint(x.geometry.values.tolist()))
# ...so have to convert back to a Geo frame after.
multipoint2 = gpd.GeoDataFrame({'name':multipoint.index}, geometry=multipoint.values, crs={'init':'epsg:4326'})
convex_hulls = gpd.GeoDataFrame(multipoint2.name, geometry=multipoint2.convex_hull, crs={'init':'epsg:4326'})
convex_hulls.to_file('/home/at15963/Dropbox/work/papers/bamber_fwf/ice_discharge/bedmachine_catchments.shp')