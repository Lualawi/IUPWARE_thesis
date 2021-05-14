#!/bin/bash

name=cesm2
model=tas_${name}_cutforelev
observation=obs_cutforelev

cdo remapcon2,../CMIP/$model.nc ../elev/elev_GMTED.nc ../elev/elev_tas.nc
cdo remapcon2,../CMIP/$observation.nc ../elev/elev_tas.nc ../elev/elev_grid_$name.nc
cdo remapcon2,../CMIP/$observation.nc ../elev/elev_GMTED.nc ../elev/elev_obs.nc

rm -f ../elev/elev_tas.nc

