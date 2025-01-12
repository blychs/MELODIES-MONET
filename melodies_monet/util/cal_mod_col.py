# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

# calculate WRF-Chem NO2 trop. columns, for further pair with satellite swath data
# contact: meng.li.atm@gmail.com
#

import numpy as np
import xarray as xr
import pandas as pd

def cal_model_columns(modobj, sp='no2'):

    """
    Calcuate model columns for each layer, to pair with satellite data
    
    Parameters
    ------
    modobj : xarray.Dataset
        model data
    sp : str
        str with the species to return
  
    Output
    ------
    modobj : xarray.Dataset
        revised model data with 'no2col' and 'localtime' added
    """

    # calculate the no2 tropospheric vertical columns and pressure from wrf-chem
    # update, mli, to be coordinated with monetio
    species    = modobj[sp]
    tb     = modobj['temperature_k']
    pb2    = modobj['pres_pa_mid'] # time,z,y,x
    dzh    = modobj['dz_m'] # time,z,y,x, the model layer thickness 
    time   = modobj.coords['time']
    modlon = modobj.coords['longitude']


    # presure: base state + PB (KSMP)
    nt, nz, ny, nx = species.shape
    #pb2  = np.zeros([nt, nz, ny, nx],dtype=np.float)
    #pb2  = pdata + pb

    # convert the perturbation potential temperature (from 300K reference) to temp
    #tb = np.zeros([nt, nz, ny, nx],dtype=np.float)
    #tb =(300.0+tdata)*((pb2/1.0e5)**0.286)


    # --- initialize arrays
    # no2 columns for each layer
    speciescol     = np.zeros([nt, nz, ny, nx], dtype = np.float32)
    # temporary array
    value      = np.zeros([nt, ny, nx], dtype = np.float32)

    # average between 13:00 and 14:00 localtime
    localtm    = np.zeros([nt,ny,nx], dtype='datetime64[s]')
    tdlist     = np.zeros([ny], dtype=np.int32)
    tdlt       = np.zeros([ny, nx], dtype='timedelta64[ms]')

    for xx in range(nx):
         tdlist[:]  = np.array(modlon[:,xx]/15.0).astype(int)
         tdlt[:,xx] = pd.TimedeltaIndex(tdlist, 'h')

    for tt in range(nt):
        localtm[tt,:,:] = pd.Timestamp(time.values[tt]) + tdlt[:,:]

    # --- calculate NO2 columns by layers
    # convert to ppm
    species = species / 1000.0
    for vl in range(nz):
        ad = pb2[:,vl,:,:] * (28.97e-3)/(8.314*tb[:,vl,:,:])
        #zh = ((ph[:,vl+1,:,:] + phb[:,vl+1,:,:]) - (ph[:,vl,:,:]+phb[:,vl,:,:]))/9.81
        value[:,:,:]= species[:,vl,:,:]*dzh[:,vl,:,:]*6.022e23/(28.97e-3)*1e-10*ad[:,:,:] # timex y x x
        speciescol[:,vl,:,:] = value[:,:,:]

    # add to model
    #modobj['PB2'] = xr.DataArray(pb2,dims=["time", "z", "y","x"]) # change from "time" to "t"
    modobj['localtime'] = xr.DataArray(localtm, dims=["t","y", "x"])
    modobj[f'{sp}col'] = xr.DataArray(speciescol,dims=["t", "z", "y","x"])

    return modobj

def cal_model_total_column(modobj, sp):
    """
    Calcuate model columns for each layer, to pair with satellite data
    
    Parameters
    ------
    modobj : xarray.Dataset
        model data
    sp : str
        str with the species to return
  
    Output
    ------
    modobj : xarray.Dataset
        revised model data with 'col' and 'localtime' added
    """

    data = cal_model_columns(modobj, sp)
    modobj[f'{sp}totalcol'] = modobj[f'{sp}col'].sum(dim='z', keep_attr=True)
    return modobj
