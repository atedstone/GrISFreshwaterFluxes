""" 

Downgrade solid ice discharge

"""

import xarray as xr
import pandas as pd
import numpy as np

folder_path = '/home/at15963/Dropbox/work/papers/bamber_fwf/outputs_Jan2018/'

dataset = xr.open_dataset(folder_path + 'FWF17.v3.nc', chunks={'TIME':6})

dataset.solid_ice.resample(TIME='1AS').sum().resample(TIME='1MS').apply(lambda x: x / 12)

new = (dataset.solid_ice.resample(TIME='1AS', keep_attrs=True).sum(dim='TIME') / 12).resample(TIME='1MS', keep_attrs=True).ffill()

# 2016 values don't ffill to the end of the year, so do this manually.
v2016 = new.sel(TIME='2016-01-01').load().values
new_new = xr.DataArray(np.repeat(np.array([v2016]), 11, axis=0), dims=['TIME', 'Y', 'X'], coords={'TIME':pd.date_range('2016-02-01','2016-12-01', freq='1MS'), 'Y':dataset.Y, 'X':dataset.X})

concat = xr.concat([new, new_new], dim='TIME')

dataset['solid_ice'] = concat
dataset.solid_ice.attrs['description'] = 'the monthly discharge data are mean annual values divided by 12.'

print('Writing dataset . . . ')
dataset.to_netcdf(folder_path + 'FWF17.v3_a.nc')