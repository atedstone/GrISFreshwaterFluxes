"""
First run solid_ice_discharge.py in workspace
"""

## For paper: calculate % of ice-sheet-wide total discharged by King outlets
stats_king_annual = sid_glaciers_monthly \
	.filter(items=enderlin_names_for_king) \
	.sum(axis=1).resample('1AS').sum()

stats_all_annual = sid_glaciers_monthly.sum(axis=1).resample('1AS').sum()

# This is the same or almost the same as stats_all_annual
stats_enderlin_annual = sid_glaciers_monthly.filter(items=enderlin_cols) \
	.sum(axis=1).resample('1AS').sum()

stats_comp_annual = (100 / stats_all_annual['2000':'2012']) * stats_king_annual['2000':'2012']
stats_comp_mean = stats_comp_annual.mean()



## Basins

# n.b. basin12 = enderlin higher than rignot, so some cancelling here
major_east = ['basin13', 'basin45', 'basin12']
# basin 13 runs from -42.81548,61.82051 to -40.44215,65.17088
# and 45 from -34.59111,66.60769 to -33.16263,67.57917

east_tot_rignot = basin_q_rignot.filter(items=major_east)['2000':'2009']
east_tot_enderlin = basin_q_enderlin.filter(items=major_east)['2000':'2009']
(east_tot_rignot.sum(axis=1) - east_tot_enderlin.sum(axis=1)).mean()


major_sw = ['basin14', 'basin15', 'basin16' ,'basin17']
# These basins correspond to south and south west coast, from -45.05625,61.52447 to -49.62371,64.83733
sw_tot_rignot = basin_q_rignot.filter(items=major_sw)['2000':'2009']
sw_tot_enderlin = basin_q_enderlin.filter(items=major_sw)['2000':'2009']
(sw_tot_rignot.sum(axis=1) - sw_tot_enderlin.sum(axis=1)).mean()


