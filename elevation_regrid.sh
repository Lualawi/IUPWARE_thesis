#!/bin/bash

name=cesm2
model=tas_${name}_cutforelev
observation=obs_cutforelev

cdo remapcon2,../CMIP/$model.nc ../elev/elev_GMTED.nc ../elev/elev_tas_$name.nc
cdo remapbil,../elev/$observation.nc ../elev/elev_tas_$name.nc ../elev/elev_grid_$name.nc


rm -f ../elev/elev_tas.nc

