"""

Combine runoff, discharge and masks into single netcdf

"""

import xarray as xr
import georaster

folder_path = '/home/at15963/Dropbox/work/papers/bamber_fwf/outputs_Jan2018/'

runoff = xr.open_dataset(folder_path + 'FWF17_runoff_RACMO2.3p2.nc')
discharge = xr.open_dataset(folder_path + 'FWF17_solidice_RACMO2.3p2.nc')

coords = {'Y':runoff.Y, 'X':runoff.X}

LSMGr = georaster.SingleBandRaster(folder_path + 'mask_LSMGr_filled.tif')
LSMGr_xr = xr.DataArray(LSMGr.r, 
	coords=coords, 
	dims=['Y', 'X'], 
	encoding={'dtype':'int8', 'zlib':True, '_FillValue':-99})
LSMGr_xr.name = 'LSMGr'
LSMGr_xr.attrs['long_name'] = 'Hole-filled Greenland land mass mask'
LSMGr_xr.attrs['units'] = 'none'
LSMGr_xr.attrs['grid_mapping'] = 'polar_stereographic'

basins = georaster.SingleBandRaster(folder_path + 'outflow_basins.tif')
meta = basins.ds.GetMetadata()
basins_xr = xr.DataArray(basins.r, 
	coords=coords, 
	dims=['Y', 'X'], 
	encoding={'dtype':'int8', 'zlib':True, '_FillValue':-99})
basins_xr.name = 'ocean_basins'
basins_xr.attrs['long_name'] = 'ID number of oceanographic basin which each coastal pixel drains into.'
basins_xr.attrs['units'] = 'none'
basins_xr.attrs['grid_mapping'] = 'polar_stereographic'
basins_xr.attrs['basins'] = meta['basins']
basins_xr.attrs['history'] = meta['history']


merged = xr.merge([runoff, discharge, LSMGr_xr, basins_xr])
merged.attrs = runoff.attrs
merged.attrs.pop('history')
merged.attrs['title'] = 'Monthly freshwater fluxes to the ocean across the Arctic, 1958-2016'

merged.to_netcdf(folder_path + 'FWF17.v3.nc')
