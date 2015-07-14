# -*- coding: utf-8 -*-
"""Read OSKAR binary files from python."""

import sys
import struct
import collections
import numpy


class OskarBinary(object):

    """Class providing an interface to OSKAR binary data files.

    see:
        http://www.oerc.ox.ac.uk/~ska/oskar2/OSKAR-Binary-File-Format.pdf

    TODO:
        - Split data reading from indexing to be able to deal with very large
          files and make reading a sub-set of the data faster.
    """

    class DataType:
        Char, Int, Single, Double, _, Complex, Matrix, _ = range(8)

    class Group:
        _, Standard, _, Settings, RunInfo, _, _,\
            Sky, _, Spline, Element, VisHeader, VisBlock = range(13)

    class Standard:
        _, DateTime, Version, UserName, WorkingDir = range(5)

    class Settings(object):
        Path = 1
        File = 2

    class RunInfo(object):
        Log = 1

    def __init__(self, file_name):
        """Constructor."""
        self.file_name = file_name
        self.file_handle = open(file_name)
        self.bin_ver = 0
        self.data = collections.OrderedDict()
        self.read()

    def __del__(self):
        """Destructor."""
        self.file_handle.close()

    def read_header(self):
        """Read header."""
        f = self.file_handle
        name = f.read(9)
        if name[0:8] != 'OSKARBIN':
            raise RuntimeError('Not a valid OSKAR binary file.')
        bin_ver = struct.unpack('B', f.read(1))[0]
        if not (bin_ver == 1 or bin_ver == 2):
            raise RuntimeError('The class can only read OSKAR binary '
                               'format version 1 or 2.')
        self.bin_ver = bin_ver

        # Version 1: header information.
        if bin_ver == 1:
            endian = struct.unpack('B', f.read(1))[0]
            svoid = struct.unpack('B', f.read(1))[0]
            sint = struct.unpack('B', f.read(1))[0]
            slong = struct.unpack('B', f.read(1))[0]
            sfloat = struct.unpack('B', f.read(1))[0]
            sdouble = struct.unpack('B', f.read(1))[0]
            patch = struct.unpack('B', f.read(1))[0]
            minor = struct.unpack('B', f.read(1))[0]
            major = struct.unpack('B', f.read(1))[0]
            other = struct.unpack('B', f.read(1))[0]
        # Version 2: read remaining reserved space.
        else:
            _ = f.read(64 - 10)

    @staticmethod
    def is_set(x, n):
        """Checks if a flag is set (value of bit n in byte x)."""
        return x & 2**n != 0

    def read_block_header(self, block_index):
        """."""
        f = self.file_handle

        element_size = struct.unpack('B', f.read(1))[0]
        chunk_flags = struct.unpack('B', f.read(1))[0]
        data_type = struct.unpack('B', f.read(1))[0]
        group = struct.unpack('B', f.read(1))[0]
        tag = struct.unpack('B', f.read(1))[0]
        index = struct.unpack('i', f.read(4))[0]
        block_size = struct.unpack('l', f.read(8))[0]

        if group not in self.data:
            self.data[group] = collections.OrderedDict()
        if tag not in self.data[group]:
            self.data[group][tag] = collections.OrderedDict()
        if index not in self.data[group][tag]:
            self.data[group][tag][index] = collections.OrderedDict()

        block = self.data[group][tag][index]
        block['group'] = block_size
        block['tag'] = block_size
        block['index'] = block_size
        block['number'] = block_index
        block['element_size'] = element_size
        block['chunk_flags'] = chunk_flags
        block['flag_endian'] = self.is_set(chunk_flags, 5)
        block['flag_crc'] = self.is_set(chunk_flags, 6)
        block['flag_extended'] = self.is_set(chunk_flags, 7)
        block['data_type'] = data_type
        block['block_size'] = block_size

        return block

    def read_block_data(self, block):
        """."""
        f = self.file_handle

        # Data size of the block payload.
        data_size = block['block_size']
        if block['flag_crc']:
            data_size -= 4

        # Read the block payload.
        if self.is_set(block['data_type'], self.DataType.Char):
            name = 'char'
            n = data_size
            data = f.read(data_size)

        elif self.is_set(block['data_type'], self.DataType.Int):
            name = 'int'
            n = data_size / block['element_size']
            data = struct.unpack('i' * n, f.read(data_size))

        elif self.is_set(block['data_type'], self.DataType.Single):
            if self.is_set(block['data_type'], self.DataType.Matrix):
                if self.is_set(block['data_type'], self.DataType.Complex):
                    name = 'single complex matrix'
                    n = data_size / block['element_size'] * 2 * 4
                else:
                    name = 'single matrix'
                    n = data_size / block['element_size'] * 4
            else:
                if self.is_set(block['data_type'], self.DataType.Complex):
                    name = 'single complex'
                    n = data_size / block['element_size'] * 2
                else:
                    name = 'single'
                    n = data_size / block['element_size']
            data = struct.unpack('f' * n, f.read(data_size))

        elif self.is_set(block['data_type'], self.DataType.Double):
            if self.is_set(block['data_type'], self.DataType.Matrix):
                if self.is_set(block['data_type'], self.DataType.Complex):
                    name = 'double complex matrix'
                    n = data_size / block['element_size'] * 2 * 4
                else:
                    name = 'double matrix'
                    n = data_size / block['element_size'] * 4
            else:
                if self.is_set(block['data_type'], self.DataType.Complex):
                    name = 'double complex'
                    n = data_size / block['element_size'] * 2
                else:
                    name = 'double'
                    n = data_size / block['element_size']
            data = struct.unpack('d ' * n, f.read(data_size))

        else:
            raise RuntimeError('ERROR: Unknown binary data type detected.')

        # Add the data block into the block dictionary.
        block['data_type_name'] = name
        block['data_length'] = n
        block['data'] = numpy.squeeze(data)

        if block['flag_crc']:
            # TODO(BM) implement CRC check. e.g. http://goo.gl/IfyyOO
            f.read(4)

    def read_data(self):
        """."""
        f = self.file_handle
        block_id = 0
        while f.read(3) == b'TBG':
            block = self.read_block_header(block_id)
            self.read_block_data(block)
            block_id += 1

    def read(self):
        """."""
        self.read_header()
        self.read_data()

    def date_time(self):
        gid = self.Group.Standard
        tid = self.Standard.DateTime
        if gid in self.data and tid in self.data[gid]:
            assert len(self.data[gid][tid]) == 1, \
                'Expecting only one standard group, date-time tag!'
            return self.data[gid][tid][0]['data']

    def user(self):
        gid = self.Group.Standard
        tid = self.Standard.UserName
        if gid in self.data and tid in self.data[gid]:
            assert len(self.data[gid][tid]) == 1, \
                'Expecting only one standard group, user tag!'
            return self.data[gid][tid][0]['data']

    def settings(self):
        gid = self.Group.Settings
        tid = self.Settings.File
        if gid in self.data and tid in self.data[gid]:
            assert len(self.data[gid][tid]) == 1, \
                'Expecting only one standard group, settings tag!'
            return self.data[gid][tid][0]['data']

    def print_summary(self):
        for group_id in self.data:
            group_data = self.data[group_id]
            for tag_id in group_data:
                tag_data = group_data[tag_id]
                for index in tag_data:
                    block = tag_data[index]
                    print '[%03i]' % block['number'],
                    block_id = '%i.%i.%i' % (group_id, tag_id, index)
                    print '%-9s' % block_id,
                    if block['flag_crc']:
                        print 'crc',
                    print ''


class OskarVis(OskarBinary):

    """."""

    class VisHeader:
        TelescopePath = 1
        NumVisBlockTags = 2
        FlagAutoCorrelation = 3
        FlagCrossCorrelation = 4
        VisDataType = 5
        CoordDataType = 6
        MaxTimes = 7
        NumTimes = 8
        MaxChannels = 9
        NumChannels = 10
        NumStations = 11
        PolarisationType = 12
        PhaseCentreCoordType = 21
        PhaseCentre = 22
        StartFrequency = 23
        FrequencyIncrement = 24
        ChannelBandwidth = 25
        StartTime = 26
        TimeInterval = 27
        TimeIntegration = 28
        TelescopeLon = 29
        TelescopeLat = 30
        TelescopeAlt = 31
        StationX = 32
        StationY = 33
        StationZ = 34

    class VisBlock:
        Dims = 1
        AutoCorrelation = 2
        CrossCorrelation = 3
        UU = 4
        VV = 5
        WW = 6

    def __init__(self, file_name):
        OskarBinary.__init__(self, file_name)
        # super(OskarVis, self).print_summary()
        if not self.bin_ver == 2:
            raise ValueError("Only OSKAR binary format version-2.0 files "
                             "can be read by this class.")

        # Make local copies of visibility header variables.
        vis_header = self.data[self.Group.VisHeader]
        assert len(vis_header) == 26, \
            'Expecting the visibility header to have 26 tags!'
        self.block_length = vis_header[self.VisHeader.MaxTimes][0]['data']
        self.num_times = vis_header[self.VisHeader.NumTimes][0]['data']
        self.num_channels = vis_header[self.VisHeader.NumChannels][0]['data']
        self.num_stations = vis_header[self.VisHeader.NumStations][0]['data']
        self.num_baselines = self.num_stations * (self.num_stations - 1) / 2
        self.num_blocks = int(numpy.ceil(float(self.num_times) /
                                         self.block_length))
        #
        # block_dims = self.data[self.Group.VisBlock][self.VisBlock.Dims]
        # for index in block_dims:
        #     print index, block_dims[index]['data']

    def uvw(self, flatten=False):
        # FIXME(BM) handle channels?
        group = self.Group.VisBlock
        tag_uu = self.VisBlock.UU
        tag_vv = self.VisBlock.VV
        tag_ww = self.VisBlock.WW
        uu = numpy.empty((self.num_baselines, self.num_times), dtype='f8')
        vv = numpy.empty((self.num_baselines, self.num_times), dtype='f8')
        ww = numpy.empty((self.num_baselines, self.num_times), dtype='f8')
        for index in range(0, self.num_blocks):
            block_dims = self.data[group][self.VisBlock.Dims][index]['data']
            block_times = block_dims[2]
            block_time_start = block_dims[0]
            block_baselines = block_dims[4]
            assert block_baselines == self.num_baselines, \
                "Data dimension mismatch"
            assert block_times <= self.block_length, \
                "Invalid block length ?!."
            uu_block = self.data[group][tag_uu][index]['data']
            uu_block = uu_block[0:block_baselines * block_times]
            uu_block = uu_block.reshape((block_baselines, block_times))
            uu[:, block_time_start:block_time_start + block_times] = uu_block
            vv_block = self.data[group][tag_vv][index]['data']
            vv_block = vv_block[0:block_baselines * block_times]
            vv_block = vv_block.reshape((block_baselines, block_times))
            vv[:, block_time_start:block_time_start + block_times] = vv_block
            ww_block = self.data[group][tag_ww][index]['data']
            ww_block = ww_block[0:block_baselines * block_times]
            ww_block = ww_block.reshape((block_baselines, block_times))
            ww[:, block_time_start:block_time_start + block_times] = ww_block
        if flatten:
            uu = uu.flatten()
            vv = vv.flatten()
            ww = ww.flatten()
        return uu, vv, ww


    def print_summary(self):
        print 'No. times     : %i' % self.num_times
        print 'No. channels  : %i' % self.num_channels
        print 'No. baselines : %i' % self.num_baselines
        return
        for group_id in self.data:
            group_data = self.data[group_id]
            for tag_id in group_data:
                tag_data = group_data[tag_id]
                for index in tag_data:
                    block = tag_data[index]
                    print '[%03i]' % block['number'],
                    block_id = '%i.%i.%i' % (group_id, tag_id, index)
                    print '%-9s' % block_id,
                    group_name = ''
                    if group_id == self.Group.VisHeader:
                        group_name = 'VisHeader'
                    if group_id == self.Group.VisBlock:
                        group_name = 'VisBlock'
                    print '%-15s' % group_name,
                    if block['flag_crc']:
                        print 'crc',
                    print ''


if __name__ == "__main__":

    if not len(sys.argv) == 2:
        print 'Usage: $ python read_oskar_vis.py <filename>'
        sys.exit(1)

    file_name = sys.argv[1]

    vis = OskarVis(file_name)
    # vis.print_summary()
    # print vis.date_time()
    # print vis.user()
    # print vis.settings()
    # vis.print_summary()
    # #print vis.num_blocks()

    uu, vv, ww = vis.uvw(flatten=True)
    import matplotlib.pyplot as plt
    uu = numpy.concatenate((uu, -uu))
    vv = numpy.concatenate((vv, -vv))
    plt.plot(uu, vv, '.')
    plt.show()


