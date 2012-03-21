#!/usr/bin/env python
'''
Converts XSQ files to FASTQ format
'''

import os
import sys
import gzip
import multiprocessing

from xsqutils import XSQFile

try:
    from eta import ETA
except:
    ETA = None


def pretty_number(n):
    count_l = list(str(n))
    for i in range(len(count_l))[::-3][1:]:
        count_l.insert(i + 1, ',')
    return ''.join(count_l)


def xsq_list(filename, count=False, minreads=-1):
    xsq = XSQFile(filename)
    print 'Tags: '
    for tag in xsq.tags:
        t = xsq.tags[tag]
        if t.is_colorspace:
            print '    %s[cs/%s]' % (tag, t.prefix)
        else:
            print '    %s[nt]' % (tag,)
    print ''
    print 'Samples: '

    try:
        for sample in xsq.get_samples():
            desc = xsq.get_sample_desc(sample)

            if count:
                readcount = xsq.get_read_count(sample)
                if readcount > minreads:
                    pn = pretty_number(readcount)

                    if desc:
                        print '    %s (%s) %s' % (sample, desc, pn)
                    else:
                        print '    %s %s' % (sample, pn)
            else:
                if desc:
                    print '    %s (%s)' % (sample, desc)
                else:
                    print '    %s' % (sample, )
    except KeyboardInterrupt:
        pass

    xsq.close()


def xsq_info(filename):
    xsq = XSQFile(filename)
    xsq.dump(xsq.hdf.root.RunMetadata)
    xsq.close()


def _xsq_convert_region(filename, sample, region, tags, outname):
    out = gzip.open(outname, 'w')
    xsq = XSQFile(filename)

    for name, seq, quals in xsq.fetch_region(sample, region, tags):
        if suffix:
            out.write('@%s%s\n%s\n+\n%s\n' % (name, suffix, seq, ''.join([chr(q + 33) for q in quals])))
        else:
            out.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(q + 33) for q in quals])))
    xsq.close()
    out.close()
    return region


def _dump_stream(src, dest, chunk_size=4 * 1024 * 1024):  # use 4MB chunk to read/write
    while True:
        buf = src.read(chunk_size)
        dest.write(buf)
        if len(buf) < chunk_size:
            return


class Callback(object):
    def __init__(self, total):
        self.i = 0
        self.eta = ETA(total)

    def __call__(self, result=None):
        self.i += 1
        self.eta.print_status(self.i, extra=result)

    def done(self):
        self.eta.done()


#  TODO: Make this multi-process - add job queue? Or just workers?
def xsq_convert(filename, sample=None, tags=None, suffix=None, procs=1, outname='-', tmpdir='.', noz=False):
    sys.stderr.write("Converting: %s\n" % sample)

    if procs < 1:
        procs = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(procs)

    xsq = XSQFile(filename)
    regions = []
    tmpnames = []
    for region in xsq.get_regions(sample):
        regions.append(region)
        tmpnames.append(os.path.join(tmpdir, '.tmp.%s.%s.%s.fastq.gz' % (os.path.basename(filename), sample, region)))
    xsq.close()

    if ETA:
        callback = Callback(len(regions))
    else:
        callback = None

    for region, tmpname in zip(regions, tmpnames):
        pool.apply_async(_xsq_convert_region, (filename, sample, region, tags, tmpname), callback=callback)

    pool.close()
    try:
        pool.join()
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    if callback:
        callback.done()

    sys.stderr.write("Merging temp files...\n")
    if ETA:
        callback = Callback(len(regions))
    else:
        callback = None

    tmpname = os.path.join(tmpdir, '.tmp.%s.%s.%s' % (os.path.basename(outname), sample, os.getpid()))

    if outname == '-':
        out = sys.stdout
    elif noz:
        out = open(tmpname, 'w')
    else:
        out = gzip.open(tmpname, 'w')

    for tmpname in tmpnames:
        src = gzip.open(tmpname)
        _dump_stream(src, out)
        src.close()
        os.unlink(tmpname)
        if callback:
            callback()

    if out != sys.stdout:
        out.close()
        os.rename(tmpname, outname)

    if callback:
        callback.done()


def xsq_convert_all(filename, tags=None, force=False, suffix=None, noz=False, usedesc=False, minreads=0, fsuffix=None, unclassified=False, procs=1):
    xsq = XSQFile(filename)

    samples = []

    for sample in xsq.get_samples():
        fname = sample
        if not fsuffix:
            fsuffix = ''

        if usedesc:
            fname = xsq.get_sample_desc(sample)
            if not fname:
                fname = sample

        if fname == sample:
            sys.stderr.write('Sample: %s... ' % fname)
        else:
            sys.stderr.write('Sample: (%s) %s... ' % (sample, fname))

        if noz:
            outname = '%s%s.fastq' % (fname, fsuffix)
        else:
            outname = '%s%s.fastq.gz' % (fname, fsuffix)

        if force or not os.path.exists(outname):
            if sample == 'Unclassified' and not unclassified:
                sys.stderr.write(' Skipping unclassified\n')
                continue

            count = xsq.get_read_count(sample)
            if count < minreads:
                sys.stderr.write(' Too few reads (%s)\n' % count)
                continue

            samples.append((sample, outname))
        sys.stderr.write('\n')

    xsq.close()

    for sample, outname in samples:
        xsq_convert(filename, sample, tags, suffix, procs=procs, outname=outname, noz=noz)


def usage():

    print '''Usage: xsq cmd {opts} filename.xsq

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
'''
    sys.exit(1)


if __name__ == '__main__':
    cmd = None
    sample_name = None
    tags = []
    all = False
    fnames = []
    last = None
    force = False
    count = False
    suffix = None
    noz = False
    minreads = 0
    procs = 1
    usedesc = False
    fsuf = None
    unclassified = False

    for arg in sys.argv[1:]:
        if not cmd and arg in ['list', 'convert', 'info']:
            cmd = arg
        elif last == '-n':
            sample_name = arg
            last = None
        elif last == '-s':
            suffix = arg
            last = None
        elif last == '-t':
            tags .append(arg)
            last = None
        elif last == '-min':
            minreads = int(arg)
            last = None
        elif last == '-procs':
            procs = int(arg)
            last = None
        elif last == '-fsuf':
            fsuf = arg
            last = None
        elif arg in ['-t', '-n', '-s', '-min', '-fsuf', '-procs']:
            last = arg
        elif arg == '-noz':
            noz = True
        elif arg == '-c':
            count = True
        elif arg == '-f':
            force = True
        elif arg == '-a':
            all = True
        elif arg == '-desc':
            usedesc = True
        elif arg == '-unclassified':
            unclassified = True
        elif os.path.exists(arg):
            fnames.append(arg)
        else:
            print 'Unknown argument: %s' % arg

    if not cmd or not fnames:
        usage()

    for fname in fnames:
        sys.stderr.write('[%s]\n' % fname)
        if cmd == 'list':
            xsq_list(fname, count, minreads)
        elif cmd == 'info':
            xsq_info(fname)
        elif cmd == 'convert':
            if all:
                xsq_convert_all(fname, tags, force, suffix, noz, usedesc, minreads, fsuf, unclassified, procs)
            elif sample_name:
                if len(fnames) > 1:
                    sys.stderr.write('Too many files given! Must only convert one file at a time in this mode!\n\n')
                    usage()
                xsq_convert(fname, sample_name, tags, suffix, procs)
            else:
                sys.stderr.write('Missing argument! Must specify "-a" or "-n sample"\n\n')
                usage()
