#!/usr/bin/env python
'''
Converts XSQ files to FASTQ format
'''

import os
import sys
import gzip
import collections

import h5py

Tag = collections.namedtuple('Tag', 'tag is_colorspace prefix')

try:
    from eta import ETA
except:
    ETA = None


class XSQFile(object):
    def __init__(self, fname):
        self.fname = fname
        self.fileobj = h5py.File(fname, 'r')
        self._samples = []

        for sample in self.fileobj:
            if sample not in  ['RunMetadata', 'Indexing']:
                self._samples.append(sample)

        self.tags = {}
        for tag in self.fileobj['RunMetadata']['TagDetails']:
            t = self.fileobj['RunMetadata']['TagDetails'][tag]
            self.tags[tag] = Tag(tag, t.attrs['IsColorPresent'][0] == 1, t.attrs['TagSequence'][0])

    def close(self):
        self.fileobj.close()

    def get_samples(self):
        return self._samples

    def get_read_count(self, sample):
        if not sample in self.fileobj:
            raise "Invalid sample name: %s" % sample
        count = 0
        for region in self.fileobj[sample]:
            count += len(self.fileobj[sample][region]['Fragments']['yxLocation'])
        return count

    def fetch(self, sample, tags=None):
        if not tags:
            tags = [t for t in  self.tags]

        if not sample in self.fileobj:
            raise "Invalid sample name: %s" % sample

        if ETA:
            count = 0
            for region in self.fileobj[sample]:
                count += len(self.fileobj[sample][region]['Fragments']['yxLocation'])
            eta = ETA(count)
        else:
            eta = None

        # reading in one location at a time should be slower than
        # just accessing the entire table at once, but this way, we
        # avoid the memory required and this code is much simpler for
        # writing out the tags interlaced. (I also suspect that h5py reads
        # the whole table into memory anyway)

        n = 0
        for i, region in enumerate(self.fileobj[sample]):
            locations = []
            for y, x in self.fileobj[sample][region]['Fragments']['yxLocation']:
                locations.append((y, x))

            for i, (y, x) in enumerate(locations):
                if eta:
                    n += 1
                    eta.print_status(n)

                for tag in tags:
                    if self.tags[tag].is_colorspace:
                        k = 'ColorCallQV'
                    else:
                        k = 'BaseCallQV'

                    basequals = self.fileobj[sample][region][tag][k][i]

                    calls = []
                    if self.tags[tag].prefix:
                        calls.append(self.tags[tag].prefix)

                    quals = []
                    for byte in basequals:
                        call = str(byte & 0x03)
                        qual = byte >> 2

                        if qual == 63:
                            call = '.'
                            qual = 0

                        calls.append(call)
                        quals.append(qual)

                    name = '%s_%s_%s' % (int(region), y, x)
                    if len(tags) > 1:
                        name = name + ' %s' % (tag)

                    yield(name, ''.join(calls), quals)
        if eta:
            eta.done()

    def dump(self, key, parent=None, indent=0):
        if parent is None:
            parent = self.fileobj

        spaces = '  ' * indent

        print '%s[%s]' % (spaces, key)

        if 'CLASS' in parent[key].attrs and parent[key].attrs['CLASS'] == 'TABLE':
            headers = []
            maxsize = []
            for name, val in parent[key].attrs.iteritems():
                if name[:6] == 'FIELD_' and name[-5:] == '_NAME':
                    headers.append(val)
                    maxsize.append(len(val))

            values = []
            for child in parent[key]:
                values.append(child)

            for valset in values:
                for i, val in enumerate(valset):
                    if len(str(val)) > maxsize[i]:
                        maxsize[i] = len(str(val))

            sys.stdout.write(spaces + '  ')
            sys.stdout.write(' | ')
            for header, size in zip(headers, maxsize):
                sys.stdout.write(header.ljust(size))
                sys.stdout.write(' | ')
            sys.stdout.write('\n')
            for valset in values:
                sys.stdout.write(spaces + '  ')
                sys.stdout.write(' | ')
                for val, size in zip(valset, maxsize):
                    sys.stdout.write(str(val).ljust(size))
                    sys.stdout.write(' | ')
                sys.stdout.write('\n')

        else:
            if parent[key].attrs:
                for name, val in parent[key].attrs.iteritems():
                    print '%s  %s: %s' % (spaces, name, val)

            for child in parent[key]:
                self.dump(child, parent[key], indent + 1)


def xsq_list(filename):
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
        count_l = list(str(xsq.get_read_count(sample)))
        for i in range(len(count_l))[::-3][1:]:
            count_l.insert(i + 1, ',')

        print '    %s (%s)' % (sample, ''.join(count_l))
    xsq.close()


def xsq_info(filename):
    xsq = XSQFile(filename)
    xsq.dump('RunMetadata')
    xsq.close()


def xsq_convert(filename, sample=None, tags=None):
    xsq = XSQFile(filename)
    for name, seq, quals in xsq.fetch(sample, tags):
        sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(q + 33) for q in quals])))
    xsq.close()


def xsq_convert_all(filename, tags=None, force=False):
    xsq = XSQFile(filename)
    for sample in xsq.get_samples():
        sys.stderr.write('Sample: %s... ' % sample)

        outname = os.path.join(os.path.dirname(filename), '%s.fastq.gz' % sample)

        if force or not os.path.exists(outname):
            sys.stderr.write('\n')
            out = gzip.open(outname, 'w')
            for name, seq, quals in xsq.fetch(sample, tags):
                out.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(q + 33) for q in quals])))
            out.close()
        else:
            sys.stderr.write('File exists! Not overwriting without -f\n')
    xsq.close()


def usage():

    print '''Usage: xsq cmd {opts} filename.xsq

Commands:
    info      - Lists all of the data associated with the XSQ file
    list      - Lists the samples and tags (R3/F3/etc) present in the file
    convert   - Converts XSQ samples and fragments to FASTQ format
        Options:
          -a        Convert all samples (saves to sample_name.fastq.gz)
          -f        Overwrite existing files
          -n name   Convert only sample "name" (writes to stdout)
                    (can be only one)
          -t tag    Convert only this tag (can be more than one)
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

    for arg in sys.argv[1:]:
        if not cmd and arg in ['list', 'convert', 'info']:
            cmd = arg
        elif last == '-n':
            sample_name = arg
            last = None
        elif last == '-t':
            tags .append(arg)
            last = None
        elif arg in ['-t', '-n']:
            last = arg
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
        xsq_list(fname)
    elif cmd == 'info':
        xsq_info(fname)
    elif cmd == 'convert':
        if all:
            xsq_convert_all(fname, tags, force)
        elif sample_name:
            xsq_convert(fname, sample_name, tags)
        else:
            sys.stderr.write('Missing argument! Must specify "-a" or "-n sample"\n\n')
            usage()
