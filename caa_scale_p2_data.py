"""

Scale CAA 2016 values

RACMO2.3p1 @1km CAA are available for 1958-2015 only.
We need to extend the time series to end of 2016 so this means using RACMO2.3p2
at 11 km resolution, which is known to not capture all runoff.
Compute a scaling factor to apply to the 2016 RACMO2.3p2 data by comparing
the model outputs for 2011-2015.
"""

import xarray as xr
import pandas as pd
import georaster 


# Calculate scaling factor between the 2.3p2 11km versus 1km 2.3p1.

# 11km-->5km pstere data, derived from RACMO 2.3p2
arctic_5km = xr.open_dataset('/scratch/process/project_RACMO/runoff_pstere_Nov2017.nc',
	chunks={'TIME':12})
times = pd.date_range('1958-01-01', '2016-12-31', freq='1MS')
arctic_5km['TIME'] = times 
ncaa_5km_mask = georaster.SingleBandRaster('/scratch/process/project_RACMO/mask_NCAA_Icemask_arctic.tif')
scaa_5km_mask = georaster.SingleBandRaster('/scratch/process/project_RACMO/mask_SCAA_Icemask_arctic.tif')

# NCAA @ 1km
ncaa_1km = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/NCAA_20180104/runoff.1958-2015.BN_RACMO2.3p1_NCAA.MM_20180104.nc',
	decode_times=False, chunks={'time':12})
ncaa_1km_masks = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/CAA_topo_icemask_lsm_lon_lat_CAA_North_NCAA.nc')
ncaa_1km_runoff = ncaa_1km.runoff.rename({'lon':'x', 'lat':'y'})
ncaa_1km_runoff = ncaa_1km_runoff.assign_coords(x=ncaa_1km_masks['x'])
ncaa_1km_runoff = ncaa_1km_runoff.assign_coords(y=ncaa_1km_masks['y'])

# SCAA @ 1km
scaa_1km = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/runoff.1958-2015.BN_RACMO2.3p1_ZGRN11_SCAA.MM.nc',
		decode_times=False, chunks={'time':10})
scaa_1km_masks = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/CAA_topo_icemask_lsm_lon_lat_CAA_South_SCAA.nc')
scaa_1km = scaa_1km.drop('lon')
scaa_1km = scaa_1km.drop('lat')
scaa_1km_runoff = scaa_1km.runoff.rename({'lon':'x', 'lat':'y'})
scaa_1km_runoff = scaa_1km_runoff.assign_coords(x=scaa_1km_masks['x'])
scaa_1km_runoff = scaa_1km_runoff.assign_coords(y=scaa_1km_masks['y'])

# Deal with time coordinate for NCAA and SCAA
times = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')
ncaa_1km_runoff['time'] = times
scaa_1km_runoff['time'] = times


## CHECK FOR FLIP UP/DOWN!!
## Calculate time series
ncaa_5km_ts = (arctic_5km.runoff \
	.where(ncaa_5km_mask.r == 1) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X', 'Y')) * (5*5) / 1e6).to_pandas()

scaa_5km_ts = (arctic_5km.runoff \
	.where(scaa_5km_mask.r == 1) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X', 'Y')) * (5*5) / 1e6).to_pandas()

ncaa_1km_ts = (ncaa_1km_runoff \
	.where(ncaa_1km_masks.Icemask == 1) \
	.sum(dim=('x','y')) \
	.resample(time='1AS').sum() / 1e6).to_pandas()

scaa_1km_ts = (scaa_1km_runoff \
	.where(scaa_1km_masks.Icemask == 1) \
	.sum(dim=('x','y')) \
	.resample(time='1AS').sum() / 1e6).to_pandas()


ncaa_scaling = 1. / ncaa_5km_ts['2011':'2015'].mean() * ncaa_1km_ts['2011':'2015'].mean()
ncaa_tot_2016 = ncaa_5km_ts['2016'] * ncaa_scaling ## This comes out at 1.52
# Plug this value into arctic_runoff_routing.py

# Graphing of SCAA shows no scaling factor needed for this area.

