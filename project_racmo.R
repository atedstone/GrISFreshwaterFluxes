## Project RACMO 2.3 data


library(raster)
system('mkdir -p /scratch/process/Rdump')
rasterOptions(tmpdir=paste0('/scratch/process/Rdump'))

new_proj = "+init=epsg:3413"

# Get lat/lon coordinates to use throughout operation
lats = raster('/scratch/L0data/RACMO/RACMO2.3_GRN11_masks.nc', 
		varname='lat')
lats_1d = getValues(lats)
lons = raster('/scratch/L0data/RACMO/RACMO2.3_GRN11_masks.nc', 
		varname='lon')
lons_1d = getValues(lons)

# Define function to fill nodata gaps
# https://gis.stackexchange.com/questions/181011/fill-the-gaps-using-nearest-neighbors
fill.na = function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
}  

## Define the projection function
project_layer = function(layer, round=FALSE){
	z = getValues(layer)
	# Create lon, lat, runoff
	xyz = cbind(lons_1d, lats_1d, z)
	colnames(xyz) = c('lon', 'lat', 'val')
	xyz = as.data.frame(xyz)
	# x and y here are the column headings contained in xyz
	coordinates(xyz) = ~lon+lat
	crs(xyz) = '+init=epsg:4326'

	# Project coordinates to new projection
	xyz_new_proj = spTransform(xyz, CRS(new_proj))

	# Create empty raster with same res as RACMO nominal resolution
	rast = raster(ext=extent(xyz_new_proj), crs=new_proj, resolution=11000)
	# Put irregular points onto raster
	rasOut = rasterize(xyz_new_proj, rast, xyz_new_proj$val)
	# Interpolate to 5km
	new_5km = raster(ext=extent(xyz_new_proj), crs=new_proj, resolution=5000)
	rasOut5k = resample(rasOut, new_5km, method='bilinear')
	# Filter/replace NaNs
	r2 = focal(rasOut5k, w=matrix(1,3,3), fun=fill.na, pad=TRUE, na.rm=FALSE)
	# Round it
	if(round == TRUE){
		rasOut = round(rasOut)
	}

	return(r2)
}


### Project the runoff data
nbands = 708
fname_store = vector('list', nbands)
for(i in 1:nbands){
	print(i)

	# Load in runoff for specified month
	runoff = raster('/scratch/L0data/RACMO/Nov2017/runoff.1958-2016.BN_RACMO2.3p2_FGRN11_11km.MM.nc', 
		varname='runoff', band=i)

	r2 = project_layer(runoff, round=FALSE)

	fname = paste('/scratch/process/project_RACMO/runoff_b', i, '.tif', sep='')
	fname_store[[i]] = fname
	writeRaster(r2, fname, overwrite=TRUE)
}

stacked = stack(fname_store)
# This doesn't give us a nice Time dimension, just integers...
# but I don't think it's possible with writeRaster
writeRaster(stacked, '/scratch/process/project_RACMO/runoff_pstere_Nov2017.nc', overwrite=TRUE, 
	format='CDF', varname='runoff', xname='X', yname='Y', zname='TIME', 
	longname='Runoff mm we')




### Masks
# layers = list('LSM_noGrIS', 'icecon', 'icemask_no_GrIS', 'Gr_land', 
#  	'GrIS_mask', 'GrIS_caps_mask', 'topography') 

# store = vector('list', 6)
# fname_store = vector('list', 6)
# for(i in 1:6){
# 	print(i)

# 	# Load in runoff for specified month
# 	mask = raster('/scratch/L0data/RACMO/Nov2017/FGRN11_masks.nc', 
# 		varname=layers[[i]])

# 	if(layers[[i]] == 'topography'){
# 		m = project_layer(mask, round=FALSE)
# 	}
# 	else{
# 		m = project_layer(mask, round=TRUE)
# 	}
# 	# turn NAs to zeros
# 	m[is.na(m)] = 0

# 	fname = paste('/scratch/process/project_RACMO/mask_', layers[[i]], '.tif', sep='')
# 	writeRaster(m, fname, overwrite=TRUE)
# }