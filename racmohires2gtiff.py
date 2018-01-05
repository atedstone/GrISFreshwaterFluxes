"""
gdal_translate NetCDF:"runoff.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc":runoffcorr -b 7 -of vrt testb7.vrt

gdalwarp -r average /scratch/L0data/RACMO/Nov2017/testb7.vrt /scratch/tmp/runoff_b7.tif


HOWEVER there is something weird about the vrt/tif output by gdal_translate, qgis doesn't like the look of it and it doesn't mosaic into the bigger tif with any values

It isn't picking up the geo-referencing ---- will have to write out new tiffs myself :(

MAKE SURE TO RUN using the conda environment: geospatial, NOT pygeoprocess

---

Could set the target extent of the vrt manually:
gdal_translate NetCDF:"runoff.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc":runoffcorr -b 7 -a_srs "+init=epsg:3413" testb7_2.tif

-a_srs "+init=epsg:3413"
-a_ullr ulx uly lrx lry
"""

import xarray as xr
import georaster
import pandas as pd
import gdal
import numpy as np
import subprocess
import os

PROCESS_DIR = os.environ['PROCESS_DIR']

do_gris = False
do_ncaa = False
do_scaa = False
do_mosaic = False
do_masks = False
do_grounded_ice = True

### Greenland 1 km
if do_gris:
	runoff = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/runoff.1958-2016.BN_RACMO2.4_FGRN11_GrIS.MM.nc',
		decode_times=False, chunks={'time':10})
	times = pd.date_range('1958-01-01', '2016-12-31', freq='1MS')
	runoff['time'] = times

	proj4 = '+init=epsg:3413'
	# extent and geotrans grabbed from gdalinfo
	extent = [-639956, -655596, 856044, -3355596]
	geotrans = (-639956,1000,0,-655596,0,-1000)

	n = 1
	for date in times:
		print(date)
		# Import runoff data grid from netcdf
		r_month = np.flipud(runoff.runoffcorr.sel(time=date).values.squeeze())
		georaster.simple_write_geotiff(
			PROCESS_DIR + 'GrIS_1km_runoff_month_%s.tif' % n,
			r_month,
			geotrans,
			proj4=proj4,
			dtype=gdal.GDT_Float32
			)
		n += 1

	runoff = None



### NCAA 
if do_ncaa:
	runoff = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/NCAA_20180104/runoff.1958-2015.BN_RACMO2.3p1_NCAA.MM_20180104.nc',
		decode_times=False, chunks={'time':10})
	# Load masks in order to get coordinates - not in provided runoff file!
	ncaa_masks = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/CAA_topo_icemask_lsm_lon_lat_CAA_North_NCAA.nc')
	times = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')
	#runoff['time'] = times

	proj4 = '+init=epsg:3413'
	# extent and geotrans had to be calculated as not available in ncdf
	extent = [ncaa_masks.x[0], ncaa_masks.y[-1]-1000, ncaa_masks.x[-1]-1000, ncaa_masks.y[0]]
	geotrans = (ncaa_masks.x[0],1000,0,ncaa_masks.y[-1]-1000,0,-1000)

	n = 0
	for date in times:
		print(date)
		# Import runoff data grid from netcdf
		r_month = np.flipud(runoff.runoff.isel(time=n).values.squeeze())
		georaster.simple_write_geotiff(
			PROCESS_DIR + 'NCAA_1km_runoff_month_%s.tif' % n,
			r_month,
			geotrans,
			proj4=proj4,
			dtype=gdal.GDT_Float32
			)
		n += 1

	runoff = None



### SCAA 
if do_scaa:
	runoff = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/runoff.1958-2015.BN_RACMO2.3p1_ZGRN11_SCAA.MM.nc',
		decode_times=False, chunks={'time':10})
	times = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')
	#runoff['time'] = times

	proj4 = '+init=epsg:3413'
	# extent and geotrans had to be calculated as not available in ncdf
	extent = [runoff.x[0], runoff.y[-1]-1000, runoff.x[-1]-1000, runoff.y[0]]
	geotrans = (runoff.x[0],1000,0,runoff.y[-1]-1000,0,-1000)

	n = 0
	for date in times:
		print(date)
		# Import runoff data grid from netcdf
		r_month = np.flipud(runoff.runoff.isel(time=n).values.squeeze())
		georaster.simple_write_geotiff(
			PROCESS_DIR + 'SCAA_1km_runoff_month_%s.tif' % n,
			r_month,
			geotrans,
			proj4=proj4,
			dtype=gdal.GDT_Float32
			)
		n += 1

	runoff = None



if do_mosaic:
	times = pd.date_range('1958-01-01', '2016-12-31', freq='1MS')
	n = 1
	for time in times:

		print(time)

		cmd = 'gdalwarp -r average'
		into_fn = '%s/runoff_b%s.tif' %(PROCESS_DIR, n)

		# Add GrIS
		fn = '%s/GrIS_1km_runoff_month_%s.tif' %(PROCESS_DIR, n)
		subprocess.call('%s %s %s' %(cmd, fn, into_fn), shell=True)

		if n <= 695:
			# Add NCAA
			# Use a cutline because otherwise the NCAA domain overlaps the GrIS domain
			fn = '%s/NCAA_1km_runoff_month_%s.tif' % (PROCESS_DIR, n-1)
			subprocess.call('%s -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/NCAA_bounds.shp %s %s' %(cmd, fn, into_fn), shell=True)

			# Add SCAA
			fn = '%s/SCAA_1km_runoff_month_%s.tif' % (PROCESS_DIR, n-1)
			subprocess.call('%s %s %s' %(cmd, fn, into_fn), shell=True)

		n += 1


	#gdalwarp -r average /scratch/process/project_RACMO/GrIS_1km_runoff_month_7.tif /scratch/process/project_RACMO/runoff_b7.tif
	#gdalwarp -r average -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/NCAA_bounds.shp /scratch/process/project_RACMO/NCAA_1km_runoff_month_6.tif /scratch/process/project_RACMO/runoff_b7.tif


	# Ideally need a consistency check that correct month of data is being substituted?
	# Also, are extents lining up precisely? Check for grid shifts.


"""
Masks/DEMs required:
For routing...
	DEM of ice surfaces
	ice mask
	land mask (can be combined)

For analysis...
	Need to check with JB



"""
if do_masks:

	import collections

	masks_11km = ['mask_icemask', 'mask_Geopotential']

	# Create a copy of each 11km mask to mosaic into later
	for mask in masks_11km:
		original_fn = '%s/%s.tif' %(PROCESS_DIR, mask)
		into_fn = '%s/%s_mosaic.tif' %(PROCESS_DIR, mask)
		subprocess.call('cp %s %s' %(original_fn, into_fn), shell=True)
	

	## Specify region masks to mosaic in		
	# Only Ice and Topo are available - not land masses
	# Make sure they are listed in same order as destinations in masks_11km
	region_masks = collections.OrderedDict( 
	GrIS= {'fn': 'Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc', 
			 'layers': ['Promicemask', 'Topography'] },
	NCAA= {'fn': 'CAA_topo_icemask_lsm_lon_lat_CAA_North_NCAA.nc',
			 'layers': ['Icemask', 'Topography'] },
	SCAA= {'fn': 'CAA_topo_icemask_lsm_lon_lat_CAA_South_SCAA.nc',
			 'layers': ['Icemask', 'Topography'] }
	)

	proj4 = '+init=epsg:3413'

	# Iterate region-by-region, mosaicing into domain-wide masks
	for key,vals in reversed(region_masks.items()):

		if key == 'GrIS':
			# The X and Y fields in the masks nc are broken, use runoff nc file instead
			nc = xr.open_dataset('%s/%s' %('/scratch/L0data/RACMO/Nov2017', 'runoff.1958-2016.BN_RACMO2.4_FGRN11_GrIS.MM.nc'), decode_times=False)
			extent = [float(nc.lon[0]), float(nc.lat[-1]-1000), float(nc.lon[-1]-1000), float(nc.lat[0])]
			geotrans = (float(nc.lon[0]),1000,0,float(nc.lat[-1]-1000),0,-1000)
			nc = None
			nc = xr.open_dataset('%s/%s' %('/scratch/L0data/RACMO/Nov2017', vals['fn']))
		else:
			# Open region nc
			nc = xr.open_dataset('%s/%s' %('/scratch/L0data/RACMO/Nov2017', vals['fn']))
			# extent and geotrans have to be calculated as not available in ncdf
			extent = [float(nc.x[0]), float(nc.y[-1]-1000), float(nc.x[-1]-1000), float(nc.y[0])]
			geotrans = (float(nc.x[0]),1000,0,float(nc.y[-1]-1000),0,-1000)

		print(key)
		print(extent)
		print(geotrans)

		n = 0 # counter for indexing into masks_11km
		# Iterate for each layer in region
		for layer in vals['layers']:
			# Write a geotiff so that GDAL can use it
			data = np.flipud(nc[layer].values.squeeze())
			tif_fn = PROCESS_DIR + '/mask_%s_%s.tif' %(key, layer)
			georaster.simple_write_geotiff(
				tif_fn,
				data,
				geotrans,
				proj4=proj4,
				dtype=gdal.GDT_Int16
				)

			# Now use gdalwarp to add to mosaic tif
			if layer == 'Topography':
				cmd = 'gdalwarp -r average'
			else:
				cmd = 'gdalwarp -r max'
			if key == 'NCAA':
				cmd += ' -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/NCAA_bounds.shp '
			dest_fn = '%s/%s' %(PROCESS_DIR, masks_11km[n] + '_mosaic.tif')
			full_cmd = '%s %s %s' %(cmd, tif_fn, dest_fn)
			print('**CMD: ' + full_cmd)
			subprocess.call(full_cmd, shell=True)

			# Also generate stand-alone versions with full arctic extent (helps with isolating GrIS etc)
			cmd = 'gdalwarp -overwrite -te -1780479.825 -3989808.111 1979520.175 -64808.111 -tr 5000 5000 '
			gen_fn = PROCESS_DIR + '/mask_%s_%s' %(key, layer)
			full_cmd = '%s %s.tif %s_arctic.tif' %(cmd, gen_fn, gen_fn)
			print('**CMD: ' + full_cmd)
			subprocess.call(full_cmd, shell=True)

			n += 1
		nc = None



	
## Export Greenland grounded ice mask
if do_grounded_ice:
	gris_masks_complete = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/grounded_ice_mask/Icemask_Topo_Iceclasses_lon_lat_average_1km.nc',
		decode_times=False)

	# The X and Y fields in the masks nc are broken, use runoff nc file instead
	nc = xr.open_dataset('%s/%s' %('/scratch/L0data/RACMO/Nov2017', 'runoff.1958-2016.BN_RACMO2.4_FGRN11_GrIS.MM.nc'), decode_times=False)
	extent = [float(nc.lon[0]), float(nc.lat[-1]-1000), float(nc.lon[-1]-1000), float(nc.lat[0])]
	geotrans = (float(nc.lon[0]),1000,0,float(nc.lat[-1]-1000),0,-1000)
	proj4 = '+init=epsg:3413'

	# Write a geotiff so that GDAL can use it
	data = np.flipud(gris_masks_complete.grounded_ice.values.squeeze())
	tif_fn = PROCESS_DIR + '/mask_GrIS_grounded_ice.tif' 
	georaster.simple_write_geotiff(
		tif_fn,
		data,
		geotrans,
		proj4=proj4,
		dtype=gdal.GDT_Int16
		)

	# Write to an arctic-wide geotiff so we can use it in arctic_runoff_routing.py
	cmd = 'gdalwarp -overwrite -te -1780479.825 -3989808.111 1979520.175 -64808.111 -tr 5000 5000 '
	gen_fn = PROCESS_DIR + '/mask_GrIS_grounded_ice' 
	full_cmd = '%s %s.tif %s_arctic.tif' %(cmd, gen_fn, gen_fn)
	print('**CMD: ' + full_cmd)
	subprocess.call(full_cmd, shell=True)
