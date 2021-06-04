#!/bin/bash

name=cesm2
realm=r1i1p1f1
model1=tas_${name}_${realm}
model2=tasLut_${name}_${realm}
observation=Tair_merge

cdo remapbil,../CMIP/$observation.nc ../CMIP/$model1.nc ../CMIP/${model1}_obs.nc
cdo remapbil,../CMIP/$observation.nc ../CMIP/$model2.nc ../CMIP/${model2}_obs.nc
cdo sellonlatbox,-180,180,-56,84 ../CMIP/$model1.nc ../CMIP/tas_${name}_cutforelev.nc
cdo sellonlatbox,-180,180,-90,90 ../CMIP/$model1.nc ../CMIP/tas_${name}_180.nc