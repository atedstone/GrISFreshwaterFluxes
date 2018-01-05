"""
Calculate annual GrIS and Greenland PGIC runoff sums, save to txt in order
to use for runoff-discharge interpolation in solid_ice_discharge.py.

Use unrouted runoff.

"""


import xarray as xr
import pandas as pd

ds = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/runoff.1958-2016.BN_RACMO2.4_FGRN11_GrIS.MM.nc', decode_times=False, chunks={'time':1})
times = pd.date_range('1958-01-01', '2016-12-31', freq='1MS')
ds['time'] = times
# Rename lon and lat to x and y (i.e. projected nomenclature AND matching masks file)
ds = ds.rename({'lon':'x', 'lat':'y'})

masks = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc')
# Lat/lon coordinates of masks are broken, so sub in the coordinates from the runoff dataset
masks['x'] = ds.x
masks['y'] = ds.y

# Masks grounded dataset as provided by Brice on email, late Dec 2017 
masks_grounded = xr.open_dataset('/scratch/L0data/RACMO/Nov2017/grounded_ice_mask/Icemask_Topo_Iceclasses_lon_lat_average_1km.nc', decode_times=False)
masks_grounded['x'] = ds.x
masks_grounded['y'] = ds.y


# GrIS-only (i.e. promice=3, grounded ice)
ann_sum_gris = ds.runoffcorr \
	.where(masks.Promicemask == 3) \
	.where(masks_grounded.grounded_ice == 1) \
	.sum(dim=('y', 'x')) \
	.resample('1AS', dim='time', how='sum') / 1.0e6

# GrIS and peripheral (as used in work flow)
ann_sum_perip = ds.runoffcorr
	.where(masks.Promicemask > 0) \
	.where(masks.Promicemask < 3) \
	.sum(dim=('y', 'x')) \
	.resample('1AS', dim='time', how='sum') / 1.0e6

to_save = pd.concat({'GrIS grounded':ann_sum_gris.to_pandas(), 'Periph.':ann_sum_perip.to_pandas()}, axis=1)
to_save.to_csv('~/Dropbox/work/papers/bamber_fwf/outputs_Jan2018/greenland_runoff_contributions.csv')