# -*- coding: utf-8 -*-
"""
Created on Tue May 15 12:40:44 2018

Print out the variables encoded in an OEDGE netcdf file.

@author: Jacob
"""

import netCDF4 as nc

if __name__ == '__main__':

    #ncname='C:\\Users\\Jacob\\Documents\\d-167196-ext-c-b1.nc'
    #ncname='C:\\Users\\Jacob\\Documents\\projects\\d-167196\\d-167196-ext-bg-bc1.nc'
    #ncname='C:\\Users\\Jacob\\Documents\\projects\\3dlim\\colprobe-a8.nc'
    ncname='/mnt/c/Users/Shawn/Documents/d3d_work/DIVIMP Runs/Z0-167196-001.nc'

    z=nc.Dataset(ncname,'r')
    print('file format: ',z.file_format)
    print('disk format: ',z.disk_format)
    print('')

    fmt='{:24}{:15}{:11}{:42}{}'
    print(fmt.format('TAG','UNITS','DTYPE','DIMENSIONS','DESCRIPTION'))
    print('')
    for m in z.variables.keys():
        desc=z.variables[m].long_name
        try:
            dim=[str(d) for d in z.variables[m].dimensions][0]
        except:
            dim=''
        dtype=z.variables[m].dtype.str
        if hasattr(z.variables[m], 'units'):
            unit=z.variables[m].units
        else:
            unit=''
        print(fmt.format(m, unit, dtype, dim, desc))
