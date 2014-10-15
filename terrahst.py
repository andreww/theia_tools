#!/usr/bin/env python

import os
import math
import numpy as np

def read_terra_hst(filename):
    """Read an optionally gzip compressed hst file from TERRA

    Returns numpy arrays containing data from the formatted
    hst (history) file written by TERRA. All data is parsed
    but currently only a subset is returned. Currently returned
    data is:
       nr: number of radial layers (integer)
       header: a list of strings of header information
       iter_num: a numpy array of iteration numbers corresponding to each 
                 stored step (not all time steps are stored)
       time: a numpy array of the time (in yr) of each recoreded step
       rshl: a numpy array of radii for the rn+1 radial layers
       tav: a (rshl, nsteps) array of layer avarage temperatures for
            each layer at each recorded time step 
       htop:
       hbot:
       hrad:
       heat: 
    """

    # Deal with compressed files.
    if (os.path.splitext(filename)[1] == '.gz'):
        import gzip
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'r')

    # Empty lists for data storage. Note that 
    # these are named as per the TERRA variables
    # that are written (except with _something
    # added when the variable means something to 
    # python.
    header = []
    propr = []
    rshl = []
    iter_num = []
    time = []
    tstep = []
    step = []
    fnrm = []
    rnrm = []
    divnorm_cutdiv = []
    unrm = []
    ekin = []
    htop = []
    hbot = []
    hrad = []
    heat = []
    hshr = []
    hvol = []
    hnet = []
    tnrm = []
    rat = []
    rnu = []
    tbalc_sctar = []
    tav = []
   
    # Simple state machine to read the data (as charaters)
    # into the lists. Better to do the conversion to numpy 
    # arrays in one go afterwards.
    header_num = 7
    for line in f:
        if header_num == 7:
            # First line - just nr
            nr = int(line)
            # For counting lines in the timesteps and header.
            header_block_lines = 4 + int(math.ceil((float(nr)+1.0)/6.0))
            header_block_lines_tb = header_block_lines
            running_block_lines = int(math.ceil((20.0 + float(nr) + 1.0) 
                                                 / 6.0))
            running_block_lines_tb = running_block_lines
            header_num = 6
        elif header_num >= 3:
            # Four lines of header text
            header.append(line)
            header_num = header_num - 1
        elif header_num == 2:
            # First block of data - embedded in the header
            if header_block_lines_tb > header_block_lines - 4:
                # The propr array (over 4 lines)
                propr.extend(line.split())
            elif header_block_lines_tb > 0:
                rshl.extend(line.split())
            header_block_lines_tb = header_block_lines_tb - 1
            if header_block_lines_tb == 0:
                header_num = header_num - 1
        elif header_num == 1:
            # tscale 
            tscale = float(line)
            header_num = header_num - 1
        elif header_num == 0:
            # Now into the main blocks of data
            if running_block_lines_tb == running_block_lines:
                # First line of a new block
                these_temps = []
                words = line.split()
                iter_num.append(words[0])
                time.append(words[1])
                tstep.append(words[2])
                step.append(words[3])
                fnrm.append(words[4])
                rnrm.append(words[5])
            elif running_block_lines_tb == running_block_lines - 1:
                words = line.split()
                divnorm_cutdiv.append(words[0])
                unrm.append(words[1])
                ekin.append(words[2])
                htop.append(words[3])
                hbot.append(words[4])
                hrad.append(words[5])
            elif running_block_lines_tb == running_block_lines - 2:
                words = line.split()
                heat.append(words[0])
                hshr.append(words[1])
                hvol.append(words[2])
                hnet.append(words[3])
                tnrm.append(words[4])
                rat.append(words[5])
            elif running_block_lines_tb == running_block_lines - 3:
                words = line.split()
                rnu.append(words[0])
                tbalc_sctar.append(words[1])
                these_temps.extend(words[2:])
            else:
                # Some more temps
                these_temps.extend(line.split())
            running_block_lines_tb = running_block_lines_tb - 1
            if running_block_lines_tb == 0:
                running_block_lines_tb = running_block_lines
                tav.append(these_temps)

    f.close()

    # Convert into Numpy arrays (and ignore the stuff we 
    # don't need for now)
    iter_num = np.array(iter_num).astype(np.float)
    time = np.array(time).astype(np.float)
    rshl = np.array(rshl).astype(np.float)
    tav = np.array(tav).astype(np.float).T
    htop = np.array(htop).astype(np.float)
    hbot = np.array(hbot).astype(np.float)
    hrad = np.array(hrad).astype(np.float)
    heat = np.array(heat).astype(np.float)

    return nr, header, iter_num, time, rshl, tav, htop, hbot, hrad, heat


def plot_layertemp_iter(iter_num, rshl, tav, filename=None):
    """Create a visulaisation of the time evolution of the layer temperature"""

    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(tav, extent=[iter_num.min(),iter_num.max(),
              rshl.min(),rshl.max()],
              interpolation='nearest', aspect='auto')
    ax.set_ylabel('Radius (m)')
    ax.set_xlabel('Timestep number')
    cb = fig.colorbar(im)
    cb.set_label('Temperature (K)')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

if __name__ == "__main__":
    import argparse
    import terraheat

    parser = argparse.ArgumentParser(description=
              'Generate a graphs from a TERRA hst file.')
    parser.add_argument('hstfile', help='TERRA hst file')
    parser.add_argument('-o', '--outfile', help='Ouput graph to a file')
    parser.add_argument('--heatflux', help='Generate heatflux graph',
                        action='store_true')
    args = parser.parse_args()

    nr, header, iter_num, time, rshl, tav, htop, hbot, hrad, heat = \
       read_terra_hst(args.hstfile)

    for line in header:
        print line.rstrip('\r\n').strip()

    if args.heatflux:
        terraheat.plot_heat(time, htop, hbot, hrad, heat, filename=args.outfile)
    else: 
        plot_layertemp_iter(iter_num, rshl, tav, filename=args.outfile)
    

