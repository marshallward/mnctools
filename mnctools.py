import numpy as np
import netCDF4 as nc
import os

from version import __version__

# Globals
tile_attrs = ['tile_number', 'sNx', 'sNy', 'nSx', 'nSy', 'nPx', 'nPy']


def collate(tile_fnames, output_fname, partition=None):

    # Use a sample tile to initialise the output fields
    fname = tile_fnames[0]
    tile = nc.Dataset(fname, 'r')

    # Force netCDF4 output format
    output_format = 'NETCDF4'

    # Create output file using tile's format
    output_nc = nc.Dataset(output_fname, 'w', format=output_format)

    # Determine collated grid dimensions
    tiled_dimsize = {'X': tile.Nx,
                     'Y': tile.Ny,
                     'Xp1': 1 + tile.Nx,
                     'Yp1': 1 + tile.Ny}

    # Copy global attributes
    output_attrs = [attr for attr in tile.ncattrs()
                    if not attr in tile_attrs]
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
    untiled_vars = {}
    tiled_vars = {}
    buffered_vars = {}
    for v in tile.variables:
        var = tile.variables[v]
        v_out = output_nc.createVariable(v, var.dtype, var.dimensions)

        # Copy attributes
        for attr in var.ncattrs():
            v_out.setncattr(attr, var.getncattr(attr))

        # Sort tiled variables and copy untiled variables
        if any([d in var.dimensions for d in tiled_dimsize]):
            if 'T' in var.dimensions:
                buffered_vars[v] = v_out
            else:
                tiled_vars[v] = v_out
        else:
            untiled_vars[v] = v_out

    # Before closing the tile, transfer any untiled variables
    for v in untiled_vars:
        output_nc.variables[v][:] = tile.variables[v][:]

    tile.close()

    #---
    # Copy unbuffered tiled variables
    transfer_tiles(output_nc, tile_fnames, tiled_vars)

    #---
    # Process buffered variables along time axis

    # Index partitions
    # NOTE: Assumes 'T' (i.e. time) is the top axis to be partitioned
    if 'T' in output_nc.dimensions and buffered_vars:
        t_len = len(output_nc.dimensions['T'])

        # Estimate number of partitions based on available memory
        # NOTE: np.array.nbytes requires allocation, do not use
        # NOTE: To be honest, the partition forecaster isn't very good...
        if not partition:
            partition = t_len
        elif partition == 'auto':
            v_itemsize = max([output_nc.variables[v].dtype.itemsize
                              for v in buffered_vars])
            v_size = max([output_nc.variables[v].size for v in buffered_vars])
            pbs_vmem = int(os.environ['PBS_VMEM'])

            # Memory model: 80MB + array allocation
            model_vmem = (80 * 2**20) + (v_itemsize * v_size)

            # Pad memory
            partition = 1 + int(1.25 * model_vmem) / pbs_vmem

        # Determine index bounds for partitions
        t_bounds = [(i * t_len / partition, (i+1) * t_len / partition)
                    for i in range(partition)]

    else:
        t_bounds = []

    # Begin buffered tile transfer
    for ts, te in t_bounds:
        transfer_tiles(output_nc, tile_fnames, buffered_vars, ts, te)

    output_nc.close()


def transfer_tiles(output_nc, tile_fnames, tiled_vars, ts=0, te=-1):

    for v in tiled_vars:
        v_out = output_nc.variables[v]
        dims = v_out.dimensions

        # Allocate variable field
        is_buffered = False
        v_shape = list(v_out.shape)
        for i, d in enumerate(v_out.dimensions):
            if output_nc.dimensions[d].isunlimited():
                v_shape[i] = te - ts
                is_buffered = True
        field = np.empty(v_shape, dtype=v_out.dtype)

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

            # If necessary, pull out a time-buffered sample
            if is_buffered:
                # NOTE: Assumes that 'T' is the first axis
                assert var.dimensions[0] == 'T'
                var = var[ts:te, ...]

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
        if is_buffered:
            output_nc.variables[v][ts:te, ...] = field[:]
        else:
            output_nc.variables[v][:] = field[:]
