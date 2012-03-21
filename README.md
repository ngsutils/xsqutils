xsqutils
===

Converts XSQ format files to FASTQ format files.

Requires the [pytables](http://pytables.org/) library and HDF5-devel libraries to be installed.
(Note: pytables also requires numpy, numexpr, and cython)

HDF5 libraries can be downloaded from [http://www.hdfgroup.org/HDF5/](http://www.hdfgroup.org/HDF5/). 
They can also be found in the EPEL yum repository. Pytables requires HDF5 1.6.10 or better.

The installation script `install.sh` will take care of setting up a virtualenv environment and installing Pytables and it's dependencies. The virtualenv will be setup in the 'env' directory. The default `xsq` driver script expects this setup and will automatically setup the virtualenv accordingly.

When converting a sample or entire file, multiple processors can be used to speed the conversion.

Usage
---
    xsq cmd {opts} filename.xsq

    Commands:
        info      - Lists all of the data associated with the XSQ file
        list      - Lists the samples and tags (R3/F3/etc) present in the file
            Options:
              -c           Show the number of reads present for each tag
              -min {val}   Hide samples that have less than {val} reads

        convert   - Converts XSQ samples and fragments to FASTQ format
            Options:
              -a           Convert all samples (saves to sample_name.fastq.gz)
                  [-a additional options]
                  -desc          Use descriptions for the sample name
                  -f             Overwrite existing files
                  -min {val}     Skip samples that have less than {val} reads
                  -noz           Don't compress the output FASTQ files with gzip
                  -fsuf {val}    Add suffix to file name
                  -unclassified  Export "Unclassified" library (usually skipped)

              -n name        Convert only sample "name" (writes to stdout)
                             (can be only one, written uncompressed)
              -procs {val}   Use {val} number of threads (CPUs) to convert one
                             region at a time. (default 1)
              -s suffix      Append a suffix to all read names
              -t tag         Convert only this tag (can be more than one)
                             If more than one tag is given, the sequences for
                             each read will be written out together.
                               For example:
                               @read1 F3
                               F3 seq
                               +
                               F3 qual
                               @read1 R5
                               R5 seq
                               +
                               R5 qual
                               @read2
                               ...

            The default is to convert all samples and all fragments/tags.