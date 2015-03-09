
#
# run in CASA with execfile('calculate_phase_centre_direction.py')
#

#
# casapy --nogui --nologger --log2term -c calculate_phase_centre_direction.py
#

def sph2cart(lon, lat):
    import numpy as np
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return (x,y,z)

def plot_pointings_ra_dec(pointing_file):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from itertools import product, combinations
    import csv

    data = np.loadtxt(pointing_file, dtype=np.str, comments='#', delimiter=',')
    npoints = len(data)
    #npoints = 3

    ra  = np.zeros((npoints,),dtype=np.double)
    dec = np.zeros((npoints,),dtype=np.double)
    for i in range(0, npoints):

        values = data[i]
        ra[i]  = float(values[0])
        dec[i] = float(values[1])
        time_str = values[2]

    (px,py,pz) = sph2cart(ra*(np.pi/180.0), dec*(np.pi/180.0))

    fig = plt.figure(666)
    ax = fig.add_subplot(111,projection='3d')
    ax.set_aspect("equal")

    # Draw cube
    # r = [-1, 1]
    # for s, e in combinations(np.array(list(product(r,r,r))), 2):
    #     if np.sum(np.abs(s-e)) == r[1]-r[0]:
    #         ax.plot3D(*zip(s,e), color="b")

    # Draw sphere
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="r", alpha=0.3)

    # Plot data
    s = ax.scatter(px, py, pz, c='k', s=80, lw=0)
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Ra,Dec')
    ax.view_init(azim=60, elev=-140)
    #ax.view_init(azim=60, elev=-100)
    plt.savefig('pointings.png', dpi=150)
    plt.show(block=False)

def plot_pointing_az_el(pointing_file):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from itertools import product, combinations

    data = np.loadtxt(pointing_file, comments='#', delimiter=',')
    az = data[:,1]
    el = data[:,2]

    fh = open(pointing_file)
    line = fh.readline()
    fh.close()
    ra = float(line.split()[1].split(':')[1])
    dec = float(line.split()[2].split(':')[1])

    pointing_idx = int(pointing_file[-4-2:-4])

    (px,py,pz) = sph2cart(az*(np.pi/180.0), el*(np.pi/180.0))


    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    # Draw cube
    r = [-1, 1]
    for s, e in combinations(np.array(list(product(r,r,r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s,e), color="b")

    # Draw sphere
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="r", alpha=0.3)

    # Plot data
    ntimes=len(az)
    s = ax.scatter(px, py, pz, c='b', s=5, lw=0)
    i0 = ntimes/2
    s = ax.scatter(px[i0], py[i0], pz[i0], c='k', s=50, lw=0)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.set_title('Pointing %i (RA:%.3f, Dec:%.3f)' %  (pointing_idx, ra, dec), fontsize=8)
    ax.view_init(azim=az[i0]-20, elev=el[i0]-20)
    plt.savefig('pointing_%02i.png' % pointing_idx, dpi=150)
    plt.show(block=True)
    plt.close()

def generate_pointings_az_el(pointing_file, lon0, lat0, obs_length_days, ntimes):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from itertools import product, combinations

    data = np.loadtxt(pointing_file, dtype=np.str, comments='#', delimiter=',')
    npoints = len(data)

    for p in range(0, npoints):
        pointingFile = os.path.join('data', 'pointing_%02i.txt' % p)

        if os.path.exists(pointingFile): os.remove(pointingFile)

        f_handle = file(pointingFile, 'a')

        values = data[p]

        ra  = float(values[0])
        dec = float(values[1])
        time_str = values[2]
        az0  = float(values[3])
        el0 = float(values[4])
        e0 = me.epoch('UTC', time_str)
        mjd0 = e0['m0']['value']

        # Write header
        f_handle.write('# RA:%f Dec:%f\n' % (ra, dec))
        f_handle.write('# lon:%f lat:%f\n' % (lon0, lat0))
        f_handle.write('# centre time:%s %f\n' % (time_str, mjd0))
        f_handle.write('# obs length days:%f\n' % (obs_length_days))
        f_handle.write('# no. times: %i\n' % (ntimes))

        # Convert RA, Dec to Az, el for range of times.
        d0 = me.direction('J2000', qa.quantity(ra, 'deg'), qa.quantity(dec, 'deg'))
        p0 = me.position('WGS84', qa.quantity(lon0, 'deg'), qa.quantity(lat0, 'deg'),
        qa.quantity(alt, 'm'))
        mjd_start = mjd0 - obs_length_days/2.0
        mjd_inc = obs_length_days/(ntimes-1)
        az = np.zeros((ntimes,),dtype=np.double)
        el = np.zeros((ntimes,),dtype=np.double)
        for i in range(0, ntimes):
            obs_time_mjd = mjd_start + i*mjd_inc# + mjd_inc/2.0
            qt = qa.quantity(obs_time_mjd, unitname='d')
            obs_time = me.epoch('UTC', qa.time(qt, form="ymd")[0])
            time = obs_time['m0']['value']
            me.doframe(obs_time)
            me.doframe(p0)
            d = me.measure(d0, 'AZEL')
            az[i] = d['m0']['value']*(180.0/np.pi)
            el[i] = d['m1']['value']*(180.0/np.pi)
            f_handle.write('%f, %f, %f\n' % (time, az[i], el[i]));

        f_handle.close()

if __name__ == '__main__':

    import shutil
    from pprint import pprint
    import numpy as np
    import numpy.random as rand
    import os

    # ------------------------------------------------------------------
    lon0     = 21.442909
    lat0     = -30.739475
    alt      = 0.0
    t0_year  = 2015
    t0_month = 3
    t0_day   = 5
    t0_hour  = 0
    t0_min   = 0
    t0_sec   = 0
    t0_utc   = '%04i/%02i/%02i/%02i:%02i:%05.2f' % \
        (t0_year, t0_month, t0_day, t0_hour, t0_min, t0_sec)
    min_el = 60
    npointings = 2
    ntimesobs  = 5
    obs_length_days = 6.0/24.0
    outfile = os.path.join('data', 'pointings.txt')
    # -----------------------------------------------------------------

    if os.path.isdir(os.path.dirname(outfile)):
        shutil.rmtree(os.path.dirname(outfile))

    os.makedirs(os.path.dirname(outfile))
    f_handle = file(outfile, 'a')

    # Seed the random number generator
    seed = 1
    rand.seed(seed)

    t_obs_inc = obs_length_days / float(ntimesobs) # hours

    # Lon,lat of the telescope
    p0 = me.position('WGS84', qa.quantity(lon0, 'deg'), qa.quantity(lat0, 'deg'),
        qa.quantity(alt, 'm'))

    for p in range(0, npointings):
        # Genreate a random observation time for the centre of the observation.
        t_hour = rand.randint(0,24)
        t_min  = rand.randint(0,60)
        t_utc  = '%04i/%02i/%02i/%02i:%02i:%05.2f' % \
            (t0_year, t0_month, t0_day, t_hour, t_min, t0_sec)

        # Genreate an Az,El of the phase centre at the centre time.
        az = rand.randint(0,2) # Choose az of the mid-point to be either 0 or 180
        if az == 1: az = 180
        # Choose elevation of the mid-point from range (min_el < el < 90)
        el = rand.randint(min_el,90)
        # Convert az, el and time to celestial coordiantes.
        d0 = me.direction('AZEL', qa.quantity(az, 'deg'), qa.quantity(el, 'deg'))
        e = me.epoch('UTC', t_utc)
        t_mjd_utc = e['m0']['value']
        qt = qa.quantity(t_mjd_utc, unitname='d')
        qa.time(qt, form="ymd")
        me.doframe(e)
        me.doframe(p0)
        d = me.measure(d0, 'J2000')
        ra  = d['m0']['value']*(180./np.pi)
        dec = d['m1']['value']*(180./np.pi)
        stime = qa.time(qt, form=["ymd"])[0]
        if p == 0:
            f_handle.write('# RA, Dec, CASA date-time, MJD, az, el\n')
            f_handle.write('# Times, az, and el are for the mid point of the observation.\n')
            f_handle.write('')
        f_handle.write('% 15.10f, % 15.10f, %20s, %.6f, % 9.4f, % 9.4f\n' % \
            (ra, dec, stime, t_mjd_utc, az, el))

    f_handle.close()

    # for each pointing make an az el file showing the observation track.
    generate_pointings_az_el(outfile, lon0, lat0, obs_length_days, ntimesobs)
    #plot_pointings_ra_dec(outfile)
