xsqutils
===

Converts XSQ format files to FASTQ format files.

Requires the h5py library and hdf5-devel packages to be installed. It is suggested that this be done using virtualenv.

Virtualenv should be setup in the 'env' directory. HDF5 libraries should either be installed on the system or in the env/hdf5 directory. The default xsq driver script expects this setup and will automatically setup the virtualenv environment and setup the LD_LIBRARY_PATH (Linux) and DYLD_LIBRARY_PATH (Mac) environmental variables.