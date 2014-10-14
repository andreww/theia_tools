#!/usr/bin/env python

import numpy as np

def read_terra_heat(file):
  """Read an optionally compressed heat file from TERRA

     Returns five numpy arrays of summary data from a 
     TERRA heat file. The arrays contain the time (years), 
     power exiting via the surface (W), entering via 
     bottom heating (W), added via the decay of radioisotopes
     (W) and the total power (W), respectivly. 
  """

  # We can use the numpy loadtxt to read this file format, which is
  # nice - and allows file to be a file object or file name, and the 
  # file can be .gz or .bz2 compressed.
  time, htop, hbot, hrad, heat = np.loadtxt(file, unpack=True)
  return time, htop, hbot, hrad, heat

def plot_heat(time, htop, hbot, hrad, heat, filename=None):
    """Create a graph of data from a TERRA heat file

       Creates a line graph of the contributions to the power
       of the system using data that can be read from a 
       TERRA heat file. Data must be presented as numpy 
       arrays (or a format that can be converted into a numpy
       array). Setting the optional filename argument will result
       in the graph being written to a file.
    """

    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, htop, 'g-', label='htop')
    ax.plot(time, hbot, 'b-', label='hbot')
    ax.plot(time, hrad, 'r-', label='hrad')
    ax.plot(time, heat, 'k-', label='heat')
    ax.legend()
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel('Heat (W)')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=
              'Generate a graph from a TERRA heat file.')
    parser.add_argument('histfile', help='TERRA heat file')
    parser.add_argument('-o', '--outfile', help='Ouput graph to a file')
    args = parser.parse_args()

    time, htop, hbot, hrad, heat = read_terra_heat(args.histfile)
    plot_heat(time, htop, hbot, hrad, heat, filename=args.outfile)
