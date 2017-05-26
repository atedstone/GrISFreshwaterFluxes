"""

Check/verify routed runoff

"""

import xarray as xr

masks = xr.open_dataset('/scratch/process/RACMO2.3_GRN11_masks_pstere.nc')
unrouted = xr.open_dataset('/scratch/process/RACMO2.3_GRN11_runoff_monthly_1958-2015_pstere.nc')
routed = xr.open_dataset('/home/at15963/Dropbox/routed_1958_1963.nc')
routed = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/outputs/FWF17_runoff.nc', chunks={'TIME':12})

# Check GrIS-only ice runoff over 5-year mean corresponding to 2012 paper
routed.runoff_ice \
	.where(masks.Greenland_all == 1) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) \
	.load() \
	.rolling(TIME=5) \
	.mean()

routed.runoff_ice \
	.where(masks.Greenland_all == 1) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) \
	.load() \

	.mean()

routed.runoff_tundra \
	.where(masks.Greenland_all) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) \
	.rolling(TIME=5) \
	.mean()


# Need to also compared unrouted and routed fluxes


# Check ice runoff from GrIS
r5y_un = (unrouted.runoff \
	.sel(TIME=slice('1958','1964')) \
	.where(masks.GrIS_mask) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) * (5*5) / 1.0e6) \
	.rolling(TIME=5).mean()

routed.runoff_ice \
	.sel(TIME=slice('1958','1963')) \
	.where(masks.Greenland_all) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) 


# Check whole domain
unrouted_sum = unrouted.runoff \
	.sel(TIME=slice('1958','1963')) \
	.where(masks.all_ice) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) * (5*5) / 1.0e6

routed_sum = routed.runoff_ice \
	.sel(TIME=slice('1958','1963')) \
	.resample('1AS', dim='TIME', how='sum') \
	.sum(dim=('X','Y')) 

"""
As of 9 Feb - these values don't quite match, routed values are consistently 3-7 km3 more than unrouted.

NOT during the routing from ice margin to coast. Checked this in arctic_runoff_routing.py, 
runoff quantities in ice-margin==1 and then coastal grid are the same at c. 4dp.

Looks like the routing toolbox is not preserving mass.
Alternatively, when I'm identifying the outlet points by using the ice sheet margin, 
I could somehow be including additional non-outlet points.
The toolbox doesn't seem to have a way of identifying outlet points. 
HOWEVER - the difference is only c. 1% and routed systematically higher than 
unrouted - so won't have too much of an impact on long-term time series? Could 
check this I suppose.

"""


## Do a check with original unprojected, unrouted data
masks_geo = xr.open_dataset('/scratch/L0data/RACMO/RACMO2.3_GRN11_masks.nc')
runoff_geo = xr.open_dataset('/scratch/L0data/RACMO/RACMO2.3_GRN11_runoff_monthly_1958-2015.nc',
	decode_times=False, chunks={'time':12})
# Apply time dimension in format understandable by xarray
times = pd.date_range('1958-01-01', '2015-12-01', freq='1MS')
runoff_geo['time'] = times

r1ygeo = ((runoff_geo.runoff * masks_geo.gridarea.values / 1.0e6) \
	.where(masks_geo.GrIS_mask.values == 1) \
	.resample('1AS', dim='time', how='sum') \
	.sum(dim=('lon','lat')) ) \
	.load() 

brunoff = xr.open_dataset('/scratch/process/RACMO2.3_GRN11_runoff_monthly_1958-2015_bamber_TEST.nc')
grid_proj = pyproj.Proj('+proj=sterea +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
lon = runoff_geo.LON.values.flatten()
lat = runoff_geo.LAT.values.flatten()
latlon = np.vstack((lat, lon)).T
(x, y) = grid_proj(lon, lat)

for x1, y1, lon1, lat1 in zip(x, y, lon, lat):
	geo = runoff_geo.runoff.sel(lon=lon1, lat=lat1, method='nearest').values()
	proj = brunoff.runoff.sel(x=x1, y=y1, method='nearest').values()
	print('%s : %s' %(geo, proj))