#!/usr/bin/env python

def is_set(x, n):
    return x & 2**n != 0

# TODO oskarpy: function to read (group.tag.index)
# TODO oskarpy vis function to read named quantities such as
#   - uvw, amp, flag ...
#   - vis summary function.

def tag_name(group, tag):
    # Standard group
    if group == 1:
        root = 'general.'
        if tag == 1:
            return root+'file_date'
        elif tag == 2:
            return root+'oskar_version'
        elif tag == 3:
            return root+'username'
        elif tag == 4:
            return root+'working_directory'
        return root

    # Settings group
    if group == 3:
        root = 'settings.'
        if tag == 1:
            return root+'path'
        elif tag == 2:
            return root+'values'
        return root

    # Run info. group
    if group == 4:
        root = 'run_info.'
        if tag == 1:
            return root+'log'
        return root

    # Visibilities (deprecated)
    if group == 5:
        root = 'vis.'
        return root

    if group == 11:
        root = 'vis_header.'
        if tag == 1:
            return root+'telescope_path'
        if tag == 2:
            return root+'write_autocorrelations'
        if tag == 3:
            return root+'amp_type'
        if tag == 4:
            return root+'times_per_block'
        if tag == 5:
            return root+'times_total'
        if tag == 6:
            return root+'num_channels'
        if tag == 7:
            return root+'num_stations'
        if tag == 8:
            return root+'freq_start_hz'
        if tag == 9:
            return root+'freq_inc_hz'
        if tag == 10:
            return root+'channel_bandwidth_hz'
        if tag == 11:
            return root+'time_start_mjd_utc'
        if tag == 12:
            return root+'time_inc_s'
        if tag == 13:
            return root+'time_average_s'
        if tag == 14:
            return root+'phase_centre'
        if tag == 15:
            return root+'telescope_centre'
        if tag == 16:
            return root+'station_x_offset_ecef'
        if tag == 17:
            return root+'station_y_offset_ecef'
        if tag == 18:
            return root+'station_z_offset_ecef'
        return root


    if group == 12:
        root = 'vis_block.'
        if tag == 1:
            return root+'dim_start_and_size'
        if tag == 2:
            return root+'freq_ref_inc_hz'
        if tag == 3:
            return root+'time_ref_inc_mjd_utc'
        if tag == 4:
            return root+'auto_correlations'
        if tag == 5:
            return root+'cross_correlations'
        if tag == 6:
            return root+'baseline_ww'
        if tag == 7:
            return root+'baseline_vv'
        if tag == 8:
            return root+'baseline_ww'
        if tag == 9:
            return root+'baseline_num_time_averages'
        if tag == 10:
            return root+'baseline_num_channel_averages'
        return root

    return ''

if __name__ == "__main__":
    import os
    import sys
    import struct
    import binascii

    if not len(sys.argv) == 2:
        print 'Usage: $ python read_oskar_vis.py <filename>'
        sys.exit(1)

    filename = sys.argv[1]
    print 'filename    : %s' % filename
    print_data = False

    #f = open(filename, 'rb')
    f = open(filename)
    try:
        # Read the file header.
        print '[HEADER]'
        name      = f.read(9)
        bin_ver   = struct.unpack('B', f.read(1))[0]
        endian    = struct.unpack('B', f.read(1))[0]
        svoid     = struct.unpack('B', f.read(1))[0]
        sint      = struct.unpack('B', f.read(1))[0]
        slong     = struct.unpack('B', f.read(1))[0]
        sfloat    = struct.unpack('B', f.read(1))[0]
        sdouble   = struct.unpack('B', f.read(1))[0]
        patch     = struct.unpack('B', f.read(1))[0]
        minor     = struct.unpack('B', f.read(1))[0]
        major     = struct.unpack('B', f.read(1))[0]
        other     = struct.unpack('B', f.read(1))[0]
        oskar_ver = "%i.%i.%i" % (major, minor, patch)
        f.read(44)
        assert(name == b'OSKARBIN\0')
        print 'name        :', name
        print 'binary ver  :', bin_ver
        if bin_ver == 1:
            print 'endian      :', endian
            print 'size void   :', svoid
            print 'size int    :', sint
            print 'size long   :', slong
            print 'size float  :', sfloat
            print 'size double :', sdouble
        print 'oskar ver   :', oskar_ver

        print '='*80
        print ''

        # Read the data block header
        blockidx = 0
        while f.read(3) == b'TBG':
            #if blockidx > 0: break
            print ''
            print '[BLOCK %04i]' % (blockidx)
            element_size = struct.unpack('B', f.read(1))[0]
            chunk_flags = struct.unpack('B', f.read(1))[0]
            big_endian = is_set(chunk_flags,5)
            crc = is_set(chunk_flags,6)
            extended = is_set(chunk_flags,7)
            data_type = struct.unpack('B', f.read(1))[0]
            group_id = struct.unpack('B', f.read(1))[0]
            tag_id = struct.unpack('B', f.read(1))[0]
            index = struct.unpack('i', f.read(4))[0]
            block_size = struct.unpack('l', f.read(8))[0]
            print 'Element size   :', element_size
            print 'Group ID       :', group_id
            print 'Tag ID         :', tag_id
            print 'Index          :', index
            print 'Chunk flags    :', chunk_flags
            print '  - Big endian :', big_endian
            print '  - CRC        :', crc
            print '  - Extended   :', extended
            print 'Data type      :', data_type
            print 'Block size     :', block_size
            tag_name_ = tag_name(group_id, tag_id)
            if len(tag_name_) > 0:
                print 'Tag name       :', tag_name_
            data_size = block_size
            if (crc): data_size -= 4
            print '-'*80
            if is_set(data_type, 0):
                n = data_size
                print 'DATA: <char>', n
                data = f.read(data_size)
            elif is_set(data_type, 1):
                n = data_size/element_size
                print 'DATA: <int>', n
                data = struct.unpack('i'*n, f.read(data_size))
            elif is_set(data_type, 2):
                if is_set(data_type, 6):
                    if is_set(data_type, 5):
                        n = data_size/element_size
                        n *= 2*4
                        print 'DATA: <matrix complex float>', n
                        data = struct.unpack('f'*n, f.read(data_size))
                    else:
                        n = data_size/element_size
                        n *= 4
                        print 'DATA: <matrix float>', n
                        data = struct.unpack('f'*n, f.read(data_size))
                else:
                    if is_set(data_type, 5):
                        n = data_size/element_size
                        n *= 2
                        print 'DATA: <complex float>', n
                        data = struct.unpack('f'*n, f.read(data_size))
                    else:
                        n = data_size/element_size
                        print 'DATA: <float>', n
                        data = struct.unpack('f'*n, f.read(data_size))
            elif is_set(data_type, 3):
                if is_set(data_type, 6):
                    if is_set(data_type, 5):
                        n = data_size/element_size
                        n *= 2*4
                        print 'DATA: <matrix complex double>', n, data_size
                        data = struct.unpack('d'*n, f.read(data_size))
                        #data = f.read(data_size)
                    else:
                        n = data_size/element_size
                        n *= 4
                        print 'DATA: <matrix double>', n
                        data = struct.unpack('d'*n, f.read(data_size))
                else:
                    if is_set(data_type, 5):
                        n = data_size/element_size
                        n *= 2
                        print 'DATA: <complex double>', n
                        data = struct.unpack('d'*n, f.read(data_size))
                    else:
                        n = data_size/element_size
                        print 'DATA: <double>', n
                        data = struct.unpack('d'*n, f.read(data_size))
            else:
                print 'ERROR: should not be here...'
                data = f.read(data_size)
            if group_id == 11 and tag_id == 2:
                print bool(data[0])
            print '-'*80

            #print data
            if (crc):
                # implement check: http://goo.gl/IfyyOO
                f.read(4)
            blockidx += 1

    finally:
        f.close()
