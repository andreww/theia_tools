#!/usr/bin/env python

import numpy as np

def visc_profile(nu_lid, z_lid, nu_j, lam, z_j, z, z_lid_slope_width=0.0):
    """Smoothed lower viscosity jump for base of the asthenosphere
    
    The upper boundary is defined by a lid nu_lid times more 
    viscous than the asthenosphere with a base at depth z_lid
    
    For the lower boundary we use the function of Sinha and Butler 
    (2007; JGR 112:B10406 http://dx.doi.org/10.1029/2006JB004850)
    which is parameterised by:
        * z_j: the depth (in km) of the viscosty jump at the "base" 
            of the asthernosphere
        * lam: a parameter giving the width of the jump (numbers 
            between 20 and 200 give good values), larger is smoother
        * nu_j: the radial viscosity factor: how much more viscous is the
            mantle below the jump, compared to the manlte above the jump.
    """
    z_j = z_j / np.max(z)
    z_norm = z / np.max(z)
    nu = ((1-nu_j)/2)*np.tanh(lam*(z_j-z_norm))+((nu_j+1)/2)
    
    for i in range(len(nu)):
        if z[i] <= (z_lid-z_lid_slope_width):
            nu[i] = nu_lid
        elif z[i] < (z_lid):
            frac = (z_lid - z[i])/z_lid_slope_width
            nu[i] = (1-frac)*nu[i] + frac*nu_lid
               
    nu[-1] = nu_j
    
    return nu

def write_viscosity_file(filename, header, z, nu):
    f = open(filename, 'w')
    f.write(header)
    for zi, nui in zip(z, nu):
        f.write("{:10.3f} {:10.3f}\n".format(zi, nui))
    f.close()

def make_viscosity_profile(nu_lid, z_lid, nu_j, lam, z_j, z_max, steps,
       z_lid_slope_width):

    z = np.linspace(0.0, z_max, steps)
    nu = visc_profile(nu_lid, z_lid, nu_j, lam, z_j, z, z_lid_slope_width)
    return (z, nu)

def plot_viscosity_profile(z, nu, filename=None):
    import matplotlib
    if filename is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(nu, z)
    ax.invert_yaxis()
    ax.set_ylabel('Depth (km)')
    ax.set_xlabel('Viscosity factor')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=
               "Generate a viscosity profile for TERRA.")
    parser.add_argument('--nu_lid', type=float, default=10.0,
               help="Viscosity of lithosphere")
    parser.add_argument('--z_lid', type=float, default=100.0,
               help="Depth of base of lithosphere")
    parser.add_argument('--z_lid_jump', type=float, default=50.0,
               help="Smooth the LAB jump over z_lid_jump km")
    parser.add_argument('--nu_j', type=float, default=10.0,
               help="Viscosity of lower mantle")
    parser.add_argument('--lam', type=float, default=10.0,
               help="Smoothing parameter for lower viscosity jump")
    parser.add_argument('--z_j', type=float, default=400.0,
               help="Depth of lower viscosity jump")
    parser.add_argument('--z_max', type=float, default=1000.0,
               help="Maximum depth to include")
    parser.add_argument('--steps', type=int, default=100,
               help="Number of depth steps")
    parser.add_argument('-o', '--outfile', help='Write ouput to a file')
    parser.add_argument('-p', '--plot', help='Create a graph',
                        action='store_true')
    args = parser.parse_args()

    z, nu = make_viscosity_profile(args.nu_lid, args.z_lid, 
                 args.nu_j, args.lam, args.z_j, args.z_max, args.steps,
                 args.z_lid_jump)

    if args.plot:
        plot_viscosity_profile(z, nu, args.outfile)
    else:
        if args.outfile is None:
            filename = "visc_fact.dat"
        else: 
            filename = args.outfile
        write_viscosity_file(filename, "TERRA viscosity input", z, nu)


