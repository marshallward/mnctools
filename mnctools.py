#!/usr/bin/env python
# coding: utf-8

import numpy as np
import netCDF4 as nc


# Globals
tiling_attrs = ['tile_number', 'sNx', 'sNy', 'nSx', 'nSy', 'nPx', 'nPy']
#tiled_dims = ['X', 'Y', 'Xp1', 'Yp1']


def collate(tile_fnames, output_fname):
    # TODO: Big ass function, break this into testable pieces
    
    output_nc = nc.Dataset(output_fname, 'w')
    
    # Process the first element to initialise the output fields
    # TODO: It would be better to construct this from some manifest file, but
    # that may be a job for MITgcm
    
    fname = tile_fnames[0]
    tile = nc.Dataset(fname, 'r')
    
    # Determine collated grid dimensions
    tiled_dimsize = {'X': tile.Nx,
                     'Y': tile.Ny,
                     'Xp1': 1 + tile.Nx,
                     'Yp1': 1 + tile.Ny}
    
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
        elif d in tiled_dimsize:
            output_nc.createDimension(d, tiled_dimsize[d])
        else:
            output_nc.createDimension(d, len(dim))
    
    # Create a variable manifest
    tiled_vars = {}
    for v in tile.variables:
        var = tile.variables[v]
        
        v_out = output_nc.createVariable(v, var.dtype, var.dimensions)
        
        # Identify tiled variables, copy nontiled variables
        if any([d in var.dimensions for d in tiled_dimsize]):
            tiled_vars[v] = v_out
        else:
            output_nc.variables[v][:] = tile.variables[v][:]
    
    tile.close()

    # Copy tiled variables
    for v in tiled_vars:
        v_out = output_nc.variables[v]
        dims = v_out.dimensions
        
        # Allocate variable field
        v_shape = list(v_out.shape)
        for i, d in enumerate(v_out.dimensions):
            if tile.dimensions[d].isunlimited():
                v_shape[i] = v_out.shape[i]
        field = np.empty(v_shape)
        
        # Copy each tile to field
        for fname in tile_fnames:
            tile = nc.Dataset(fname, 'r')
            var = tile.variables[v]
            
            # Determine bounds: xs <= x < xe, ys <= y < ye
            nt = tile.tile_number - 1
            xt, yt = nt % tile.nPx, nt / tile.nPx
            
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
                field[..., ys:ye, xs:xe] = var[:]

            elif ('X' in dims or 'Xp1' in dims):
                field[..., xs:xe] = var[:]
            
            elif ('Y' in dims or 'Yp1' in dims):
                field[..., ys:ye] = var[:]
                
            else:
                # TODO: Use an exception?
                sys.exit('Error: untiled variable')
            
            tile.close()
        
        # Save field to output
        output_nc.variables[v][:] = field[:]

    output_nc.close()
