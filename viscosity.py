#!/usr/bin/env python

import os
import re
import sys
import numpy as np

def read_out(filename):
    """Extract the viscosity information from a TERRA out file

    The layer avarage viscosity is provided in a table at the
    end of TERRA output files. This function extracts this 
    data and places it in Numpy arrays for further 
    processing.
    """

    # Deal with compressed files.
    if (os.path.splitext(filename)[1] == '.gz'):
        import gzip
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'r')

    # Storage
    depth = []
    viscmin = []
    viscmax = []
    radial_factor = []
    viscmean = []

    # Regex setup
    re_summary = re.compile(
       r"\s+TEMPERATURE AND VISCOSITY MIN MAX SUMMARY")
    re_viscos = re.compile(
       r"\s+IR  DEPTH\(KM\)  VISCMIN  VISCMAX  RADIAL FACT  VISC MEAN") 
    re_sep = re.compile( r"\s+-{20}") 

    state = 3

    for line in f: 
        if state == 3: 
            # Outside summary 
            if re_summary.match(line):
                state = 2
        elif state == 2:
            # Not in the right table
            if re_viscos.match(line):
                state = 1
        elif state == 1:
            if re_sep.match(line):
                state = 0
        elif state == 0:
            if re_sep.match(line):
                state = 3
            else:
                words = line.split()
                depth.append(words[1])
                viscmin.append(words[2])
                viscmax.append(words[3])
                radial_factor.append(words[4])
                viscmean.append(words[5])

    f.close()

    depth = np.array(depth).astype(np.float)
    viscmin = np.array(viscmin).astype(np.float)
    viscmax = np.array(viscmax).astype(np.float)
    radial_factor = np.array(radial_factor).astype(np.float)
    viscmean = np.array(viscmean).astype(np.float)

    return depth, viscmin, viscmax, radial_factor, viscmean

def plot_viscosity(depth, viscosity, filename=None):
    """Create a graph of viscosity v's depth

    this can optionally be written in a file"""

    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.invert_yaxis()
    ax.plot(viscosity, depth, 'k-')
    ax.set_xlabel('Viscosity (Pa.s)')
    ax.set_ylabel('Depth (km)')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

def viscosity_table(depth, viscosity, filename=None):

    if filename is not None: 
        f = open(filename, 'w')
    else:
        f = sys.stdout

    for d, v in zip(depth, viscosity):
        f.write("{0:6.3f} {1:6.3e}\n".format(d, v))

    if filename is not None: 
        f.close()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=
              'Extract viscosity from TERRA out file.')
    parser.add_argument('terraout', help='TERRA out file')
    parser.add_argument('-o', '--outfile', help='Write ouput to a file')
    parser.add_argument('-p', '--plot', help='Create a graph', 
                        action='store_true')
    args = parser.parse_args()

    depth, viscmin, viscmax, radial_factor, viscmean = read_out(args.terraout)

    if args.plot:
        plot_viscosity(depth, viscmean, args.outfile)
    else:
        viscosity_table( depth, viscmean, args.outfile)
        
