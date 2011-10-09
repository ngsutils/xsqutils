xsqutils
===

Converts XSQ format files to FASTQ format files.

Requires the [h5py](http://code.google.com/p/h5py/) library and HDF5 libraries to be installed.
(Note: h5py also requires numpy)

It is suggested that this be done using virtualenv. HDF5 libraries can be downloaded from [http://www.hdfgroup.org/HDF5/](http://www.hdfgroup.org/HDF5/). 
They can also be found in the EPEL yum repository. The HDF5 libraries can either be installed system-wide or in a local directory. h5py requires HDF5 1.8.3 or better.

Virtualenv should be setup in the 'env' directory. HDF5 libraries should either be installed on the system or in the env/hdf5 directory. 
The default `xsq` driver script expects this setup and will automatically setup the virtualenv environment and setup the 
`LD_LIBRARY_PATH` (Linux) and `DYLD_LIBRARY_PATH` (Mac) environmental variables.
