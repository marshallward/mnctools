#!/usr/bin/env python
# coding: utf-8

import numpy as np
import netCDF4 as nc


# Globals
tiling_attrs = ['tile_number', 'sNx', 'sNy', 'nSx', 'nSy', 'nPx', 'nPy']

def collate(tile_fnames, output_fname, partition=None):
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
    buffered_vars = {}
    for v in tile.variables:
        var = tile.variables[v]
        
        v_out = output_nc.createVariable(v, var.dtype, var.dimensions)
        
        # Sort tiled variables and copy untiled variables
        if any([d in var.dimensions for d in tiled_dimsize]):
            if 'T' in var.dimensions:
                buffered_vars[v] = v_out
            else:
                tiled_vars[v] = v_out
        else:
            output_nc.variables[v][:] = tile.variables[v][:]
   
    tile.close()
    
    #---
    # Copy unbuffered tiled variables
    
    for v in tiled_vars:
        v_out = output_nc.variables[v]
        dims = v_out.dimensions
        
        # Allocate variable field
        v_shape = list(v_out.shape)
        for i, d in enumerate(v_out.dimensions):
            # NOTE: Unnecessary if 'T' is the only unlimited variable
            if output_nc.dimensions[d].isunlimited():
                v_shape[i] = v_out.shape[0]
        field = np.empty(v_shape)
        
        # Copy each tile to field
        for fname in tile_fnames:
            tile = nc.Dataset(fname, 'r')
            var = tile.variables[v]
           
            # Determine bounds: xs <= x < xe, ys <= y < ye
            # TODO: Precalculate and store the index ranges
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
        del field
    
    #---
    # Process buffered variables along time axis
    
    # Index partitions
    # NOTE: Assumes 'T' (i.e. time) is the top axis to be partitioned
    if 'T' in output_nc.dimensions:
        t_len = len(output_nc.dimensions['T'])
         
        # TODO: Might be able to calculate optimal partition from available
        # memory
        if not partition:
            partition = 1
        elif partition == 'all':
            partition = t_len
        
        t_bounds = [(i * t_len / partition, (i+1) * t_len / partition)
                    for i in range(partition)]
    else:
        t_bounds = []
    
    for ts, te in t_bounds:
        
        # Copy tiled variables
        for v in buffered_vars:
            v_out = output_nc.variables[v]
            dims = v_out.dimensions
            
            # Allocate variable field
            v_shape = list(v_out.shape)
            for i, d in enumerate(v_out.dimensions):
                if output_nc.dimensions[d].isunlimited():
                    v_shape[i] = te - ts
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
                    field[..., ys:ye, xs:xe] = var[ts:te, ...]
                
                elif ('X' in dims or 'Xp1' in dims):
                    field[..., xs:xe] = var[ts:te, ...]
                
                elif ('Y' in dims or 'Yp1' in dims):
                    field[..., ys:ye] = var[ts:te, ...]
                
                else:
                    # TODO: Use an exception?
                    sys.exit('Error: untiled variable')
                
                tile.close()
            
            # Save field to output
            output_nc.variables[v][ts:te, ...] = field[:]
            del field
    
    output_nc.close()
