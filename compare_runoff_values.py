"""
Compare GrIS runoff with different masks in order to confirm that Noel 2017 values are possible, and to explain
why FWF runoff is much higher due to peripheral ice and grounded ice.

Andrew Tedstone, 2 Jan 2018

"""
import xarray as xr
import pandas as pd
import georaster

ds = xr.open_dataset('runoff.1958-2016.BN_RACMO2.4_FGRN11_GrIS.MM.nc', decode_times=False, chunks={'time':1})
times = pd.date_range('1958-01-01', '2016-12-31', freq='1MS')
ds['time'] = times
# Rename lon and lat to x and y (i.e. projected nomenclature AND matching masks file)
ds = ds.rename({'lon':'x', 'lat':'y'})

masks = xr.open_dataset('Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc')
# Lat/lon coordinates of masks are broken, so sub in the coordinates from the runoff dataset
masks['x'] = ds.x
masks['y'] = ds.y

# Masks grounded dataset as provided by Brice on email, late Dec 2017 
masks_grounded = xr.open_dataset('grounded_ice_mask/Icemask_Topo_Iceclasses_lon_lat_average_1km.nc', decode_times=False)
masks_grounded['x'] = ds.x
masks_grounded['y'] = ds.y


# GrIS-only (i.e. promice=3, grounded ice)
ann_sum = ds.runoffcorr.where(masks.Promicemask == 3).where(masks_grounded.grounded_ice == 1).sum(dim=('y', 'x')).resample('1AS', dim='time', how='sum') / 1.0e6

# GrIS and peripheral (as used in work flow)
ann_sum_all = ds.runoffcorr.where(masks.Promicemask > 0).sum(dim=('y', 'x')).resample('1AS', dim='time', how='sum') / 1.0e6



# FWF dataset
fwf_r = xr.open_dataset('/home/at15963/Dropbox/work/papers/bamber_fwf/outputs_Nov2017/FWF17_runoff_RACMO2.3p2.nc')
mask_LSMGr = georaster.SingleBandRaster('/home/at15963/Dropbox/work/papers/bamber_fwf/outputs_Nov2017/mask_LSMGr_filled.tif')

fwf_ann_sum = fwf_r.runoff_tundra.where(mask_LSMGr.r == 1).sum(dim=('X', 'Y')).resample('1AS', dim='TIME', how='sum')


## Plot comparison
fwf_ann_sum.plot(label='FWF@5km')
fwf_ann_sum = fwf_r.runoff_ice.where(mask_LSMGr.r == 1).sum(dim=('X', 'Y')).resample('1AS', dim='TIME', how='sum')
figure()
fwf_ann_sum.plot(label='FWF@5km')
ann_sum_all.plot(label='RACMO1km_all')
ann_sum.plot(label='RACMO1km_GrIS')
legend()


## Look at non-GrIS FWF

fwf_nongris = fwf_r.runoff_ice.where(mask_LSMGr.r == 0).sum(dim=('X', 'Y')).resample('1AS', dim='TIME', how='sum')

## non-GrIS model domain results

ncaa = xr.open_dataset('runoff.1958-2015.BN_RACMO2.3p1_ZGRN11_NCAA.MM.nc', decode_times=False, chunks={'time':2})
ncaa_masks = xr.open_dataset('CAA_topo_icemask_lsm_lon_lat_CAA_North_NCAA.nc')
ncaa = ncaa.drop('lon')
ncaa = ncaa.drop('lat')
ncaa['x'] = ncaa_masks['x']
ncaa['y'] = ncaa_masks['y']
ncaa_runoff = ncaa.runoff.rename({'lon':'x', 'lat':'y'})
ncaa_runoff = ncaa_runoff.assign_coords(x=ncaa_masks['x'])
ncaa_runoff = ncaa_runoff.assign_coords(y=ncaa_masks['y'])

ncaa_runoff['time'] = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')
#ncaa = ncaa.rename({'lon':'x', 'lat':'y'})

figure(),ncaa_runoff.where(ncaa_masks.Icemask == 1).sum(dim=('x','y')).resample(time='1AS').sum().plot()
figure(),(ncaa_runoff.sel(time=slice('2004','2015')).where(ncaa_masks.Icemask == 1).resample(dim='time', freq='1AS', how='sum')).plot(col='time', col_wrap=5, vmin=0, vmax=500)



### South CAA
scaa = xr.open_dataset('runoff.1958-2015.BN_RACMO2.3p1_ZGRN11_SCAA.MM.nc', decode_times=False, chunks={'time':2})
scaa_masks = xr.open_dataset('CAA_topo_icemask_lsm_lon_lat_CAA_South_SCAA.nc')


## Uncorrupted NCAA
ncaa = xr.open_dataset('NCAA_20180104/runoff.1958-2015.BN_RACMO2.3p1_NCAA.MM_20180104.nc', decode_times=False)
#ncaa = ncaa.drop('lon')
#ncaa = ncaa.drop('lat')
ncaa['x'] = ncaa_masks['x']
ncaa['y'] = ncaa_masks['y']
ncaa_runoff = ncaa.runoff.rename({'lon':'x', 'lat':'y'})
ncaa_runoff = ncaa_runoff.assign_coords(x=ncaa_masks['x'])
ncaa_runoff = ncaa_runoff.assign_coords(y=ncaa_masks['y'])
ncaa_runoff['time'] = pd.date_range('1958-01-01', '2015-12-31', freq='1MS')
