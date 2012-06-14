#!/usr/bin/env python
# coding: utf-8

import numpy as np
import netCDF4 as nc


# Globals
tiling_attrs = ['tile_number', 'sNx', 'sNy', 'nSx', 'nSy', 'nPx', 'nPy']
tiled_dims = ['X', 'Y', 'Xp1', 'Yp1']


def collate(tile_fnames, output_fname):
    # TODO: Big ass function, break this into testable pieces
    
    output_nc = nc.Dataset(output_fname, 'w')
    
    # Process the first element to initialise the output fields
    # TODO: It would be better to construct this from some manifest file, but
    # that may be a job for MITgcm
    
    fname = tile_fnames[0]
    
    tile = nc.Dataset(fname, 'r')
    
    # Determine collated grid dimensions
    Nx, Ny = tile.Nx, tile.Ny
    sNx, sNy = tile.sNx, tile.sNy
    nPx, nPy = tile.nPx, tile.nPy
    tiled_dimsize = dict(zip(tiled_dims, (Nx, Ny, 1 + Nx, 1 + Ny)))
    
    # Copy global attributes
    output_attrs = [attr for attr in tile.ncattrs()
                    if not attr in tiling_attrs]
    for attr in output_attrs:
        output_nc.setncattr(attr, tile.getncattr(attr))
    
    # Copy collated dimensions
    for d in tile.dimensions:
        dim = tile.dimensions[d]
        if dim.isunlimited():
            output_nc.createDimension(d, None)
        elif d in tiled_dims:
            output_nc.createDimension(d, tiled_dimsize[d])
        else:
            output_nc.createDimension(d, len(dim))
    
    # Initialise variables, copy initial values
    
    fields = {}
    for v in tile.variables:
        var = tile.variables[v]
        dims = var.dimensions
        
        v_out = output_nc.createVariable(v, var.dtype, var.dimensions)
        
        for attr in var.ncattrs():
            v_out.setncattr(attr, var.getncattr(attr))
        
        # Pre-allocate output fields
        
        # Fix unlimited dimension sizes if necessary
        v_shape = list(v_out.shape)
        for i, d in enumerate(dims):
            if tile.dimensions[d].isunlimited():
                v_shape[i] = var.shape[i]
        fields[v] = np.empty(v_shape)
    
    tile_transfer(tile, fields)
    
    tile.close()
    
    # Copy the rest of the tile contents to the global fields
    for fname in tile_fnames[1:]:
        tile = nc.Dataset(fname, 'r')

        transfer_tile(tile, fields)
        tile.close()
    
    for v in output_nc.variables:
        var = output_nc.variables[v]
        var[:] = fields[v]
    
    output_nc.close()


def transfer_tile(tile, fields):
    
    for v in tile.variables:
        var = tile.variables[v]
        dims = var.dimensions
        
        # Determine bounds: xs <= x < xe, ys <= y < ye
        if any([d in var.dimensions for d in tiled_dims]):
            
            # Tile index
            nt = tile.tile_number - 1
            xt, yt = nt % nPx, nt / nPx
            
            xs = tile.sNx * xt
            xe = xs + tile.sNx
            ys = tile.sNy * yt
            ye = ys + tile.sNy
            
            # Extend range for boundary grids
            if 'Xp1' in dims:
                xe += 1
            if 'Yp1' in dims:
                ye += 1
         
        # Transfer tile to the collated field
        if ('X' in dims or 'Xp1' in dims) and ('Y' in dims or 'Yp1' in dims):
            fields[v][..., ys:ye, xs:xe] = var[:]
        
        elif ('X' in dims or 'Xp1' in dims):
            fields[v][..., xs:xe] = var[:]
        
        elif ('Y' in dims or 'Yp1' in dims):
            fields[v][..., ys:ye] = var[:]
            
        else:
            fields[v] = tile.variables[v][:]
