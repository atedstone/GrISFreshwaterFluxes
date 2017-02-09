import georaster
import pyproj
import osr

grid_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

im = georaster.SingleBandRaster('/home/at15963/Dropbox/work/papers/bamber_fwf/mask_oc_ed.bmp')

# (top left x, w-e cell size, 0, top left y, 0, n-s cell size (-ve))
trans = (-800000., 5000., 0., -600000., 0., -5000.)

georaster.simple_write_geotiff(
	'/home/at15963/Dropbox/work/papers/bamber_fwf/grisonly_ocean_basins.tif',
	im.r, 
	trans,
	proj4=grid_proj.srs,
	dtype=georaster.gdal.GDT_UInt16
	)
