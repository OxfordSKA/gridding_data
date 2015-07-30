#!/usr/bin/env python

def get_radii(telescope, delimiter=' '):
    """
    Returns a list of station radii for the given telescope model.

    Arguments
    ---------
    telescope : string
        Path of telescope model
    delimiter : char [default = ',']
        Delimiter of layout files within the telescope model

    Return
    ------
    List of radii of each of the child stations.
    """
    import os
    from .layout import max_radius
    # Obtain a list of any child (station) folders.
    stations = [st for st in os.listdir(telescope) \
        if os.path.isdir(os.path.join(telescope, st))]
    stations.sort()
    # Evaluate the radius of the each child (station) layout.
    r = []
    for s in stations:
        station_layout = os.path.join(telescope, s, 'layout.txt')
        r.append(max_radius(station_layout))
    return r


if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import os

    telescope = 'SKA1_mid_combined.tm'
    filename  = 'SKA1_mid_combined.eps'
    use_station_sizes = False

    # FIXME probably need to be smart about delmitiers
    delimiter=' '

    if not os.path.isdir(telescope):
        raise Exception(("Invalid telescope model, specified model path "
            "doesn't exist!"))

    # Load the top level layout file
    layout_file = os.path.join(telescope, 'layout.txt')
    if not os.path.isfile(layout_file):
        raise Exception(("Invalid telescope model, unable to find associated "
            "layout.txt file"))
    layout = np.genfromtxt(layout_file, dtype=np.double, delimiter=delimiter)


    # Generate the plot
    x = layout[:,0]/1000.0
    y = layout[:,1]/1000.0
    xymax = np.max((np.max(x), np.max(y)))*1.1
    if use_station_sizes == False:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        plt.scatter(x,y, marker='+', s=20, lw=1.0, alpha=1.0)
    else:
        # Obtain a list of the radii of any child stations.
        station_r = get_radii(telescope)
        if not len(station_r) == len(x):
            raise Exception('Unable to obtain station radii.')
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        for i in range(0, len(x)):
            if len(station_r) >= i:
                r = station_r[i]/1000.0
            c = plt.Circle((x[i],y[i]), r, alpha=1.0, color='b', ec='none')
            ax.add_artist(c)
    plt.title('Telescope layout : %s' % telescope)
    plt.xlabel('East [km]')
    plt.ylabel('North [km]')
    plt.xlim(-xymax, xymax)
    plt.ylim(-xymax, xymax)
    plt.grid()
    if not filename == None:
        plt.savefig(filename)
    plt.show()
