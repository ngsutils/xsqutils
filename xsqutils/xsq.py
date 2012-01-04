#!/usr/bin/env python
'''
Converts XSQ files to FASTQ format
'''

import os
import sys
import gzip

from xsqutils import XSQFile

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
    for sample in xsq.get_samples():
        desc = xsq.get_sample_desc(sample)
        
        if count:
            readcount = xsq.get_read_count(sample)
            if readcount > minreads:
                count_l = list(str(readcount))
                for i in range(len(count_l))[::-3][1:]:
                    count_l.insert(i + 1, ',')
            
                if desc:
                    print '    %s (%s) %s' % (sample, desc, ''.join(count_l))
                else:
                    print '    %s %s' % (sample, ''.join(count_l))
        else:
            print '    %s' % (sample, )
    xsq.close()


def xsq_info(filename):
    xsq = XSQFile(filename)
    xsq.dump('RunMetadata')
    xsq.close()


def xsq_convert(filename, sample=None, tags=None, suffix=None):
    xsq = XSQFile(filename)
    for name, seq, quals in xsq.fetch(sample, tags):
        if suffix:
            sys.stdout.write('@%s%s\n%s\n+\n%s\n' % (name, suffix, seq, ''.join([chr(q + 33) for q in quals])))
        else:
            sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(q + 33) for q in quals])))
    xsq.close()


def xsq_convert_all(filename, tags=None, force=False, suffix=None, noz=False):
    xsq = XSQFile(filename)
    for sample in xsq.get_samples():
        sys.stderr.write('Sample: %s... ' % sample)

        if noz:
            outname = os.path.join(os.path.dirname(filename), '%s.fastq' % sample)
            tmpname = os.path.join(os.path.dirname(filename), '.tmp.%s.fastq' % sample)
        else:
            outname = os.path.join(os.path.dirname(filename), '%s.fastq.gz' % sample)
            tmpname = os.path.join(os.path.dirname(filename), '.tmp.%s.fastq.gz' % sample)
            

        if force or not os.path.exists(outname):
            sys.stderr.write('\n')
            if noz:
                out = open(tmpname, 'w')
            else:
                out = gzip.open(tmpname, 'w')
                
            for name, seq, quals in xsq.fetch(sample, tags):
                if suffix:
                    out.write('@%s%s\n%s\n+\n%s\n' % (name, suffix, seq, ''.join([chr(q + 33) for q in quals])))
                else:
                    out.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(q + 33) for q in quals])))
            out.close()
            os.rename(tmpname, outname)
        else:
            sys.stderr.write('File exists! Not overwriting without -f\n')
    xsq.close()


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
          -f           Overwrite existing files
          -n name      Convert only sample "name" (writes to stdout)
                       (can be only one, written uncompressed)
          -noz         Don't compress the output FASTQ files with gzip
          -s suffix    Append a suffix to all read names
          -t tag       Convert only this tag (can be more than one)
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
    fname = None
    last = None
    force = False
    count = False
    suffix = None
    noz = False
    minreads = 0

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
        elif arg in ['-t', '-n', '-s', '-min']:
            last = arg
        elif arg == '-noz':
            noz = True
        elif arg == '-c':
            count = True
        elif arg == '-f':
            force = True
        elif arg == '-a':
            all = True
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            print 'Unknown argument: %s' % arg

    if not cmd or not fname:
        usage()

    if cmd == 'list':
        xsq_list(fname,count,minreads)
    elif cmd == 'info':
        xsq_info(fname)
    elif cmd == 'convert':
        if all:
            xsq_convert_all(fname, tags, force, suffix, noz)
        elif sample_name:
            xsq_convert(fname, sample_name, tags, suffix)
        else:
            sys.stderr.write('Missing argument! Must specify "-a" or "-n sample"\n\n')
            usage()
