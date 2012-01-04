import sys
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
    
    def get_sample_desc(self,sample):
        if sample not in self._samples:
            return None
        
        descidx = 0
        
        for name,val in self.fileobj['RunMetadata']['LibraryDetails'].attrs.iteritems():
            if val == 'Description':
                break
            descidx += 1
        
        for row in self.fileobj['RunMetadata']['LibraryDetails']:
            spl = sample.split('_')
            if spl[0] == row[0]:
                return row[descidx].strip()
                    

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

        # Reads each region into memory at a time by tag. Then yields the sequences
        # in order so that the tags are interlaced.
        #
        # This is slightly faster than just reading in one at a time.

        n = 0
        for i, region in enumerate(self.fileobj[sample]):
            locations = []
            for y, x in self.fileobj[sample][region]['Fragments']['yxLocation']:
                locations.append((y, x))

            vals = {}
            for tag in tags:
                vals[tag] = []
                if self.tags[tag].is_colorspace:
                    k = 'ColorCallQV'
                else:
                    k = 'BaseCallQV'

                for (y, x), basequals in zip(locations, self.fileobj[sample][region][tag][k]):
                    if eta:
                        n += 1
                        eta.print_status(n)

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

                    vals[tag].append((name, ''.join(calls), quals))

            for i in xrange(len(locations)):
                for tag in tags:
                    yield(vals[tag][i])
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
