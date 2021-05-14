#! /usr/bin/env python

import xarray as xr
import pandas as pd
import numpy as np

laf = xr.open_dataset('../LAF/LAF_obs_apr30.nc')

maxim1 = xr.ufuncs.maximum(laf.prim,laf.crop)
maxim = xr.ufuncs.maximum(maxim1, laf.urban)

lut = xr.where(laf.prim==maxim,1, np.NaN)
lut = lut.where((lut == 1) | (laf.crop != maxim), 2)
lut = lut.where((lut <= 2) | (laf.urban != maxim), 3)

lut = lut.rename('lut')

lut.to_dataset()

lut.to_netcdf('../LAF/LAF_remap_apr30.nc')