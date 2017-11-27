"""
gdal_translate NetCDF:"runoff.1958-2016.BN_RACMO2.3p2_FGRN11_GrIS.MM.nc":runoffcorr -b 7 -of vrt testb7.vrt

gdalwarp -r average /scratch/L0data/RACMO/Nov2017/testb7.vrt /scratch/tmp/runoff_b7.tif


HOWEVER there is something weird about the vrt/tif output by gdal_translate, qgis doesn't like the look of it and it doesn't mosaic into the bigger tif with any values

It isn't picking up the geo-referencing ---- will have to write out new tiffs myself :(

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

PROCESS_DIR = '/scratch/process/project_RACMO/'

do_gris = False
do_ncaa = False
do_scaa = False
do_mosaic = True

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
	runoff = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/runoff.1958-2015.BN_RACMO2.3p1_ZGRN11_NCAA.MM.nc',
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
		into_fn = '/scratch/process/project_RACMO/runoff_b%s.tif' % n

		# Add GrIS
		fn = '/scratch/process/project_RACMO/GrIS_1km_runoff_month_%s.tif' % n
		subprocess.call('%s %s %s' %(cmd, fn, into_fn), shell=True)

		if n <= 695:
			# Add NCAA
			# Use a cutline because otherwise the NCAA domain overlaps the GrIS domain
			fn = '/scratch/process/project_RACMO/NCAA_1km_runoff_month_%s.tif' % (n-1)
			subprocess.call('%s -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/NCAA_bounds.shp %s %s' %(cmd, fn, into_fn), shell=True)

			# Add SCAA
			fn = '/scratch/process/project_RACMO/SCAA_1km_runoff_month_%s.tif' % (n-1)
			subprocess.call('%s %s %s' %(cmd, fn, into_fn), shell=True)

		n += 1


	#gdalwarp -r average /scratch/process/project_RACMO/GrIS_1km_runoff_month_7.tif /scratch/process/project_RACMO/runoff_b7.tif
	#gdalwarp -r average -cutline /home/at15963/Dropbox/work/papers/bamber_fwf/NCAA_bounds.shp /scratch/process/project_RACMO/NCAA_1km_runoff_month_6.tif /scratch/process/project_RACMO/runoff_b7.tif


	# Ideally need a consistency check that correct month of data is being substituted?