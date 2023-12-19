"""
Python CM1 NetCDF Stitcher (cm1stitcher)
Author: Scott Powell, NPS (first dot last at nps dot edu)
Version: 0.2.2
Date: 18 December 2023

NOTE: See bottom of code for USER INPUT!
NOTE: See the README for installation instructions and details.

Description: Stitch together NetCDF output written individually by processors when running CM1 with MPI.
This will produce a single file for each output time.
Converts several files with file name convention cm1out_CPUNUM_OUTPUTNUM.nc to cm1out_OUTPUTNUM.nc

Known bugs and issues!!!
NOTE: If your domain is particularly large (i.e., more than ~250,000,000 cells), 
you may need to pay special attention to memory consumption depending on the system used.
NOTE: While I have attempted to parallelize parts of the code, performance does not appear
to scale with number of CPUs used. This needs to be fixed.
BUG: Variable 'yf' seems to not write to zarr. Reason unknown.
BUG: If looping over multiple output times, the code seems to use more memory rather than re-using the 
same amount needed for the first output time.
"""

import xarray as xr
import numpy as np
from functools import partial
import multiprocessing as mp
# from sys import argv
# from dask.distributed import LocalCluster, Client
# from glob import glob

# ***********************************************************************

# def writepkl(var,savename):
#     import pickle

#     with open(savename,'wb') as h:
#         pickle.dump(var,h,pickle.HIGHEST_PROTOCOL)

# def readpkl(fname):
#     import pickle

#     with open(fname,'rb') as h:
#         data = pickle.load(h)
#     return data

# ***********************************************************************

def setuparray(example_file):

    """
    Inputs: 
    fdir (string): Directory where data is written.
    exfile (string): Path to the example already-stitched NetCDF output file is written. 
                    Ideally, this would be a file written out by CM1 with output_filetype = 2.
    NOTE: This code requires at least one already stitched NetCDF output file
    as an example in order to grab the intended dimensions of the data. 
    This could be done by simply initializing the model to write cm1out_000001.nc.
    Otherwise, the user could alter the code to hard-code the coordinate dimensions
    as-needed.

    Outputs:
    stitched_data (xarray dataset): The array that will ultimately contain the stitched CM1 data.
    float_variables (list): The list of variables to write

    """

    # Load metadata from one of the NetCDF files to understand the dimensions
    metadata = xr.open_dataset(example_file)

    # List of variable names with float data type
    float_variables = [var for var in metadata.variables if metadata[var].dtype == 'float32']

    # List of variables to remove from float_variables.
    remove = ['ztop','xh','xf','yh','yf','zh','zf']   
    #     # These are already saved as coordinates. No need to rewrite.
    #     # Choose to not code how to write ztop with dimension one.
    #     # I figure user will know what ztop is or can figure it out from zh/zf anyway.
    for i in remove: float_variables.remove(i)

    return metadata, float_variables
    # return stitched_data, float_variables

# ***********************************************************************

def nearestvalue(A,B):

    """
    Inputs:
    A (float, int, etc.): Single value that we want to find in B.
    B (array of float, int, etc.): The array in which we want to find the value that is closest to A.
     
    Outputs:
    I (int): Index of the matching value in B.
    """

    return np.argmin(np.abs(A-B))

# ***********************************************************************

def process_file(processor,coordinates,dimensions,file_output,fdir):

    """
    Inputs: 
    processor (int): CPU number. int(XXXXXX) in cm1out_XXXXXX_YYYYYY.nc
    stitched_data
    float_variables
    file_output: int(YYYYYY) in cm1out_XXXXXX_YYYYYY.nc

    Outputs:
    stitched_data

    """

    file_name = fdir+f"cm1out_{processor:06d}_{file_output:06d}.nc"
    ds = xr.open_dataset(file_name)

    # Find the corresponding xh or xf indices based on xh or xf values
    xh_indices = [nearestvalue(i,coordinates[0]) for i in ds.xh.values]
    xf_indices = [nearestvalue(i,coordinates[1]) for i in ds.xf.values]

    # Find the corresponding yh or yf indices based on yh or yf values
    yh_indices = [nearestvalue(i,coordinates[2]) for i in ds.yh.values]
    yf_indices = [nearestvalue(i,coordinates[3]) for i in ds.yf.values]

    data = {}

    # Figure out which grid we are on for this variable.
    # x
    for var in dimensions.keys():
        if 'xh' in ds[var].coords:
            x_indices = xh_indices
        elif 'xf' in ds[var].coords:
            x_indices = xf_indices
        # y
        if 'yh' in ds[var].coords:
            y_indices = yh_indices
        elif 'yf' in ds[var].coords:
            y_indices = yf_indices

        # Insert the tile data into the stitched grid at the correct location
        if len(dimensions[var]) == 4:
            data[var] = ('all','all',y_indices[0],y_indices[-1]+1,x_indices[0],x_indices[-1]+1,ds[var].values.astype(np.float32))
            # stitched_data[var].values[:, :, y_indices[0]:y_indices[-1]+1, x_indices[0]:x_indices[-1]+1] = ds[var].values.astype(np.float32)
        elif len(dimensions[var]) == 3:
            data[var] = ('all',y_indices[0],y_indices[-1]+1,x_indices[0],x_indices[-1]+1,ds[var].values.astype(np.float32))
            # stitched_data[var].values[:, y_indices[0]:y_indices[-1]+1, x_indices[0]:x_indices[-1]+1] = ds[var].values.astype(np.float32)
        elif len(dimensions[var]) == 2:
            data[var] = (y_indices[0],y_indices[-1]+1,x_indices[0],x_indices[-1]+1,ds[var].values.astype(np.float32))
            # stitched_data[var].values[y_indices[0]:y_indices[-1]+1, x_indices[0]:x_indices[-1]+1] = ds[var].values.astype(np.float32)
        elif len(dimensions[var]) == 1:
            if var == 'xh' or var == 'xf':
                data[var] = (x_indices[0],x_indices[-1]+1,ds[var].values.astype(np.float32))
                # stitched_data[var].values[x_indices[0]:x_indices[-1]+1] = ds[var].values.astype(np.float32)
            elif var == 'yh' or var == 'yf':
                data[var] = (y_indices[0],y_indices[-1]+1,ds[var].values.astype(np.float32))
                # stitched_data[var].values[y_indices[0]:y_indices[-1]+1] = ds[var].values.astype(np.float32)

    ds.close()  # Close the dataset after use

    # Instead of returning the dictionary back to the head processor, this simply writes the data to disk.
    # After all processors write to disk, the individual files will be opened and combined.
    # savedir = '/thumper/users/scott.powell/code-data/research-code/CM1/stitching/'
    # savename = savedir + 'datatemp_'+str(processor).zfill(6)+'.pkl'
    # writepkl(data,savename)
    return data

# ***********************************************************************

def stitch(stuff,shape):

    L = len(shape)
    combine = np.zeros(shape)

    for i in stuff:
        if L == 4:
            combine[:,:,i[2]:i[3],i[4]:i[5]] = i[6]
        elif L == 3:
            combine[:,i[1]:i[2],i[3]:i[4]] = i[5]
        elif L == 2:
            combine[i[0]:i[1],i[2]:i[3]] = i[4]
        elif L == 1:
            combine[i[0]:i[1]] = i[2]

    return combine

# ***********************************************************************

# def writezarr(key):

#     # Grab the data for this variable.
#     stuff = [i[key] for i in data]

#     # Get the shape of this variable.
#     shape = stitched_data[key].shape

#     # Now populate the variable in the changable dataset.
#     var_to_zarr[key] = stitched_data[key]
#     # Change the values to the current time.
#     var_to_zarr[key].values = stitch(stuff,shape)

#     # Convert var_to_zarr to dask array then write in parallel?
#     # var_to_zarr = var_to_zarr.chunk({'time':1,'zf':201,'xh':170,'yh':315})

#     # Add this to the zarr file.
#     # If the directory doesn't exist, it will be written.
#     var_to_zarr.to_zarr(outdir+f"cm1out_{output_time:06d}.zarr",mode='a')

# ***********************************************************************

def execute(output_time, num_processors, num_workers, stitched_data, float_variables, overwrite, fdir, outdir):

    # First, grab the pieces we need from stitched_data.
    # We need the x and y coordinates at center (half) points (h) and edge (full) points (f)
    coordinates = (stitched_data.xh.values,stitched_data.xf.values,stitched_data.yh.values,stitched_data.yf.values)
    # We also need the dimensions for each variable.
    dimensions = {var: stitched_data[var].dims for var in float_variables}

    data = []
    # Run the processing in parallel. This should run much faster than it does.
    func = partial(process_file,coordinates=coordinates,dimensions=dimensions,file_output=output_time,fdir=fdir)
    with mp.Pool(processes=num_workers) as pool:
        for i in pool.imap_unordered(func,range(0,num_processors),chunksize=10):
            data.append(i)

    # Make a temporary changable copy of the stitched_data xarray dataset.
    var_to_zarr = stitched_data.copy()
    # Remove all the variables from it but leave the rest.
    vars = [var for var in var_to_zarr.variables]
    for var in vars: var_to_zarr = var_to_zarr.drop_vars(var)

    # def writezarr(key):

    #     # Grab the data for this variable.
    #     stuff = [i[key] for i in data]

    #     # Get the shape of this variable.
    #     shape = stitched_data[key].shape

    #     # Now populate the variable in the changable dataset.
    #     var_to_zarr[key] = stitched_data[key]
    #     # Change the values to the current time.
    #     var_to_zarr[key].values = stitch(stuff,shape)

    #     # Convert var_to_zarr to dask array then write in parallel?
    #     # var_to_zarr = var_to_zarr.chunk({'time':1,'zf':201,'xh':170,'yh':315})

    #     # Add this to the zarr file.
    #     # If the directory doesn't exist, it will be written.
    #     var_to_zarr.to_zarr(outdir+f"cm1out_{output_time:06d}.zarr",mode='a')

    # Write to zarr files.
    # with mp.Pool(processes=num_workers) as pool:
    #     pool.imap_unordered(writezarr,float_variables,chunksize=4)
    #     pool.close()
    #     pool.join()

    # Now stitch the data together and write to disk.
    for key in float_variables:
        stuff = [i[key] for i in data]
        shape = stitched_data[key].shape

        # Now populate the variable in the changable dataset.
        var_to_zarr[key] = stitched_data[key]
        # Change the values to the current time.
        var_to_zarr[key].values = stitch(stuff,shape)

        # Add this to the zarr file.
        # If the file doesn't exist, write it for the first time.
        var_to_zarr.to_zarr(outdir+f"cm1out_{output_time:06d}.zarr",mode='a')
        var_to_zarr = var_to_zarr.drop_vars(key)

    # Replace this to actually delete the now unneeded tiled files.
    if overwrite is True:
        pass 

# ***********************************************************************

if __name__ == '__main__':

    # *************************** USER INPUT ********************************

    # Directory where data lives.
    fdir = '/thumper/users/daniel.bazemore/CM1/Terrain_Half_Wind_Full/run/'

    # Directory where an example of fully stitched output file (like cm1out_000001.nc) lives.
    exfile = '/thumper/users/daniel.bazemore/CM1/Terrain_Half_Wind_Initial/run/cm1out_000001.nc'

    # Directory where we want to write the new stitched output.
    outdir = '/thumper/metdata/model/CM1/CACTI/test/'

    # Remove old individual tiled files? Set to True if you want to remove old files after stitching them together.
    # Recommended to only set this to True if you are certain that the newly created output meets your needs.
    # You can always go back and remove the files after, but you can't get them back after you do!
    # NOTE: This functionality is not built into current version, so setting this to True will have no effect.
    overwrite = False

    # Number of processors that was used to write CM1 output. Usually, this is the number of cores used to run CM1.
    num_processors = 512

    # What is the file output number to be processed?
    first_output_time = 2   # Integer: First output time where CM1 output is written by CPU. 
    #                         # Should be YYYYYY in cm1out_XXXXXX_YYYYYY.nc with YYYYYY as integer.
    last_output_time = 20   # Same as first_output_time, but for the last output time to be processed.
    # output_time = int(argv[1])

    # How many CPUs will be working on this? 
    # NOTE: Because of memory usage, you usually will not always simply want this to be the maximum number of CPUs available on a node. 
    num_workers = 12

    # *****************************END USER INPUT****************************
    
    # Get the common elements of all of the files to be created.
    skeleton, float_variables = setuparray(exfile)

    for output_time in range(first_output_time,last_output_time+1):
        execute(output_time, num_processors, num_workers, skeleton.copy(), float_variables, overwrite, fdir, outdir)
