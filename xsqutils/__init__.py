import sys
import collections
import re

import tables

Tag = collections.namedtuple('Tag', 'tag is_colorspace prefix')


def natural_sort(ar):
    to_sort = []
    for item in ar:
        spl = re.split('(\d+)', item)
        l2 = []
        for el in spl:
            try:
                n = int(el)
            except:
                n = el
            l2.append(n)
        to_sort.append((l2, item))

    to_sort.sort()
    return [x[1] for x in to_sort]


def node_children_iter(node):
    for k in node._v_children.keys():
        yield (k, node._f_getChild(k))


def get_node_attr(node, attr):
    val = node._v_attrs[attr]
    return convert_val(val[0], val.dtype)


def convert_val(val, dtype):
    if dtype in ['|S255', 'string']:
        return val.split('\x00')[0]
    elif dtype in ['uint8', 'int8', 'uint32']:
        return val
    else:
        sys.stderr.write("Unknown dtype: *%s*\n" % dtype)
        return val


def node_attr_iter(node):
    for k in node._v_attrs._v_attrnames:
        yield (k, get_node_attr(node, k))


class XSQFile(object):
    def __init__(self, fname):
        self.fname = fname
        self.hdf = tables.openFile(fname, 'r')
        self._samples = []

        for sample, node in node_children_iter(self.hdf.root):
            if sample not in  ['RunMetadata', 'Indexing']:
                self._samples.append(sample)

        self._samples = natural_sort(self._samples)

        self.tags = {}
        for tag, node in node_children_iter(self.hdf.root.RunMetadata.TagDetails):
            self.tags[tag] = Tag(tag, get_node_attr(node, 'IsColorPresent') == 1, get_node_attr(node, 'TagSequence'))

    def close(self):
        self.hdf.close()

    def get_samples(self):
        return self._samples

    def get_sample_desc(self, sample):
        if sample not in self._samples:
            return None

        descidx = -1
        for i, name in enumerate(self.hdf.root.RunMetadata.LibraryDetails.colnames):
            if name == 'Description':
                descidx = i

        if descidx == -1:
            return None

        desctype = self.hdf.root.RunMetadata.LibraryDetails.coltypes['Description']

        spl = sample.split('_')
        for cols in self.hdf.root.RunMetadata.LibraryDetails.cols:
            if convert_val(cols[0],self.hdf.root.RunMetadata.LibraryDetails.coltypes['LibraryName']) == spl[0]:
                return convert_val(cols[descidx], desctype)

    def get_read_count(self, sample):
        if not sample in self._samples:
            raise "Invalid sample name: %s" % sample
        count = 0

        for rn in self.hdf.root._f_getChild(sample)._v_children:
            region = self.hdf.root._f_getChild(sample)._f_getChild(rn)
            count += region._f_getChild('Fragments')._f_getChild('yxLocation').shape[0]
        return count

    def get_regions(self, sample):
        ar = self.hdf.root._f_getChild(sample)._v_children.keys()
        ar.sort()
        return ar

    def fetch_region(self, sample, region_name, tags=None):
        region_name_int = int(region_name)
        region = self.hdf.root._f_getChild(sample)._f_getChild(region_name)
        if not tags:
            tags = self.tags

        locations = []
        for y, x in region._f_getChild('Fragments')._f_getChild('yxLocation'):
            locations.append((y, x))

        vals = {}
        for tag in tags:
            vals[tag] = []
            if self.tags[tag].is_colorspace:
                k = 'ColorCallQV'
                bases = '0123'
                wildcard = '.'
            else:
                k = 'BaseCallQV'
                bases = 'ACGT'
                wildcard = 'N'

            for (y, x), basequals in zip(locations, region._f_getChild(tag)._f_getChild(k)[:]):
                name = '%s_%s_%s' % (region_name_int, y, x)
                if len(tags) > 1:
                    name = name + ' %s' % (tag)

                calls = []
                if self.tags[tag].prefix:
                    calls.append(self.tags[tag].prefix)

                quals = []
                for byte in basequals:
                    call = bases[byte & 0x03]
                    qual = byte >> 2

                    if qual == 63:
                        call = wildcard
                        qual = 0

                    calls.append(call)
                    quals.append(qual)

                vals[tag].append((name, ''.join(calls), quals))
                #yield (name, ''.join(calls), quals)

        for i in xrange(len(locations)):
            for tag in tags:
                yield(vals[tag][i])

    def dump_table(self, node, indent=0):
        spaces = '  ' * indent
        headers = []
        maxsize = []
        for name in node.colnames:
            headers.append(name)
            maxsize.append(len(name))

        values = {}
        value_names = []
        for cols in node.cols:
            for i, val in enumerate(cols):
                val1 = str(convert_val(val, node.coltypes[node.colnames[i]]))
                if len(val1) > maxsize[i]:
                    maxsize[i] = len(val1)
            name = convert_val(cols[0], node.coltypes[node.colnames[0]])
            values[name] = cols
            value_names.append(name)

        value_names = natural_sort(value_names)

        sys.stdout.write(spaces + '  ')
        sys.stdout.write(' | ')
        for header, size in zip(headers, maxsize):
            sys.stdout.write(header.ljust(size))
            sys.stdout.write(' | ')
        sys.stdout.write('\n')
        for value_name in value_names:
            sys.stdout.write(spaces + '  ')
            sys.stdout.write(' | ')
            for header, val, size in zip(headers, values[value_name], maxsize):
                sys.stdout.write(str(convert_val(val, node.coltypes[header])).ljust(size))
                sys.stdout.write(' | ')
            sys.stdout.write('\n')

    def dump(self, node, indent=0):
        if node is None:
            node = self.hdf.root

        spaces = '  ' * indent

        print '%s[%s]' % (spaces, node._v_name)

        if type(node) == tables.table.Table:
            self.dump_table(node, indent)

        else:
            for name, val in node_attr_iter(node):
                print '%s  %s: %s' % (spaces, name, val)

            for name, child in node_children_iter(node):
                self.dump(child, indent + 1)
