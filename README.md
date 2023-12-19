CM1 Stitcher 

<b>Why do we need this?</b>: CM1 can be run so that each processor individually writes output to disk at regular intervals. This is sometimes desirable because of the overhead required to collect data from all processors then write to output, often significantly increasing the run time pf CM1. While allowing each processor to write output can significantly increase the speed of CM1 and prevent compute time from being wasted on I/O, the result is numerous files at a given output time that must all be placed in the correct geographic location relative to one another. 

<b>What does this code do?</b> This code accesses all available files at a given output time, combines them, and writes them to a zarr data structure. zarr is a much more efficient storage method than NetCDF. Data can be quickly read from a zarr data structure and superior compression diminishes the disk space requirement.

<b>Caveats</b>: The speed of this code versus simply running CM1 to write all output to a single file at each output time has not been tested yet. It is possible, however, to modify this code to only write subsets of variables as needed, which could significantly decrease time required.

## Installation Instructions

This code requires xarray, numpy, and zarr.

To run the code, create a new python environment and execute the following command.

$$ pip install xarray numpy zarr

This command should also install all required dependencies.

## User input parameters

The file submitpy.sh is a bash script intended to be submitted to a SLURM scheduler. Replace Lines 2-12 as needed to meet requirements for your cluster. 

In stitchCM1.py, Lines 230-258 are user input parameters. These parameters are

fdir:               Directory where data lives.
exfile:             Directory where an example of fully stitched output file (like cm1out_000001.nc) lives.
outdir:             Directory where we want to write the new stitched output.
overwrite:          Overwrites existing zarr files, if they exist, for the output time being processed.
num_processors:     This is the number of CPUs used to write CM1 output, usually the same number of CPUs used to run CM1.
first_output_time:  Integer number for the first output time to be processed.
last_output_time:   Self explanatory.
num_workers:        Number of CPUs used to combine data in this code.

The code requires that there is at least one fully stitched CM1 output file in order to work. We usually let CM1 initialize with the same namelist, except for changing timax to 1 an changing output_filetype to 2 (it would be set to 3 if output was written in parallel) so that the model only writes the initial output file, then does nothing else. 

## Running the code

On a cluster with the SLURM scheduler, simply edit submitpy.sh to fit your needs and run 

$$ sbatch submitpy.sh

This code is currently only designed to run on a single node.

## Opening zarr data

Opening zarr data in python is easy with xarray! Suppose you have a file tree called cm1out_000010.zarr. The following snippet of code would load the uinterp variable from the zarr file tree:

import xarray as xr
zid = xr.read_zarr('cm1out_000010.zarr')
uinterp = zid['uinterp']