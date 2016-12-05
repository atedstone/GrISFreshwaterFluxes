import xarray as xr
import pyproj
from scipy import interpolate
import pandas as pd
import numpy as np

runoff = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_runoff_monthly_1958-2015.nc',
		decode_times=False, chunks={'time':15})
# Apply time dimension in format understandable by xarray
times = pd.date_range('1958-01-01', '2015-12-01', freq='1MS')
runoff['time'] = times

# Polar stereographic - same as GIMP DEM
pstere = pyproj.Proj('+init=EPSG:3413')
# Extract and project coordinates to grid
lon = runoff.LON.values.flatten()
lat = runoff.LAT.values.flatten()
(x, y) = pstere(lon, lat)
xy = np.vstack((x, y)).T

# Create Pstere grid at 11 km resolution
yi, xi = np.mgrid[np.min(y):np.max(y):11000, np.min(x):np.max(x):11000]

# Set up blank store
coords = {'TIME': times, 'Y': yi[:,0], 'X': xi[0,:]}
store = np.zeros((len(times), yi.shape[0], yi.shape[1]))

nt = 0
for t in times:

	print(t.strftime('%Y-%m'))

	r = runoff.runoff.sel(time=t).values
	r = r.flatten()

	# Interpolate point data onto grid
	zi = interpolate.griddata(xy, 
			r, (xi, yi), method='cubic')

	# Scaling factor - conservation of mass basic-style
	sf = np.sum(r) / np.nansum(zi)
	# Apply it
	zi_sf = zi * sf

	store[nt, :, :] = zi_sf
	nt += 1

da = xr.DataArray(store, coords=coords, dims=['TIME', 'Y', 'X'])
da.name = runoff.runoff.long_name
da.attrs['long_name'] = runoff.runoff.long_name
da.attrs['units'] = runoff.runoff.units
da.attrs['standard_name'] = runoff.runoff.standard_name

ds = xr.Dataset(da)
ds.attrs['history'] = 'Generated using atedstone/fwf/reproject_racmo.py on RACMO2.3_GRN11_runoff_monthly_1958-2015.nc provided by JLB/BN/MvdB'
ds.attrs['author'] = 'Andrew Tedstone a.j.tedstone@bristol.ac.uk'
ds.attrs['gridding method'] = 'cubic spline scipy.interpolate.griddata'
ds.to_netcdf('/home/at15963/Dropbox/work/papers/bamber_fwf/RACMO/RACMO2.3_GRN11_EPSG3413_runoff_monthly_1958-2015.nc')