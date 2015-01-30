#!/usr/bin/env python

import numpy as np

import tex2elas
import CijUtil

def read_texture_file(filename):
    """Read data from a theia texture file"""
    
    # Deal with compressed files.
    import os
    if (os.path.splitext(filename)[1] == '.gz'):
        import gzip
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'r')

    # Stuff everything into a dict and a list
    # for now. Sort this out later (we will probably 
    # want to have objects at some point
    header_data = {}
    particles = []

    header_lines = 5
    particle_header_lines = 9
    
    for line in f:
        if header_lines == 5:
            header_data['theia_lun'] = int(line)
            header_lines = header_lines - 1
        elif header_lines == 4:
            header_data['npartsallo'] = int(line)
            header_lines = header_lines - 1
        elif header_lines == 3:
            header_data['npartsused'] = int(line)
            header_lines = header_lines - 1
        elif header_lines == 2:
            header_data['n_expected_particles'] = int(line)
            header_lines = header_lines - 1
        elif header_lines == 1:
            header_data['nseen_particles'] = int(line)
            header_lines = header_lines - 1
        elif header_lines == 0:
            if particle_header_lines == 9:
                this_particle = {}
                this_particle['process_id'] = int(line)
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 8:
                this_particle['particle_id'] = int(line)
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 7:
                this_particle['old_particle_id'] = int(line)
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 6:
                this_particle['old_process_id'] = int(line)
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 5:
                this_particle['particle_class'] = line.strip()
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 4:
                this_particle['particle_position'] = np.array(
                    [line[0:12], line[12:24], line[24:36]])
                particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 3:
                this_particle['idata_count'] = int(line)
                if this_particle['idata_count'] > 0:
                    particle_header_lines = particle_header_lines - 1
                else:
                    particle_header_lines = particle_header_lines - 2
            elif particle_header_lines == 2:
                    this_particle['particle_idata'] = np.array(
                        [line.rstrip('\r\n')[i:i+12] for i in xrange(0, len(line.rstrip('\r\n')), 12)]
                        )
                    particle_header_lines = particle_header_lines - 1
            elif particle_header_lines == 1:
                this_particle['rdata_count'] = int(line)
                if this_particle['rdata_count'] > 0:
                    particle_header_lines = particle_header_lines - 1
                else:
                    particles.append(this_particle)
                    particle_header_lines = 9
            elif particle_header_lines == 0:
                    this_particle['particle_rdata'] = np.array(
                        [line.rstrip('\r\n')[i:i+14] for i in xrange(0, len(line.rstrip('\r\n')), 14)]
                        )
                    particles.append(this_particle)
                    particle_header_lines = 9
    f.close()

    return header_data, particles

class drex_particle(object):

    def __init__(self, position, pclass, rdata, idata):
        self.position = position
        self.pclas = pclass
        self.num_grains = idata[0]
        self._unpack_drex_rdata(rdata, self.num_grains)

    def olivine_cij(self, scheme='Voigt'):

        # Single crystal olivine elasticity from DRex - should
        # check that this is up to date.
        ol_cij_single = np.zeros((6,6))
        ol_cij_single[0,0] = 320.71
        ol_cij_single[1,1] = 197.25
        ol_cij_single[2,2] = 234.32
        ol_cij_single[0,1] = 69.84
        ol_cij_single[1,0] = ol_cij_single[0,1]
        ol_cij_single[2,0] = 71.22
        ol_cij_single[0,2] = ol_cij_single[2,0]
        ol_cij_single[1,2] = 74.80
        ol_cij_single[2,1] = ol_cij_single[1,2]
        ol_cij_single[3,3] = 63.77
        ol_cij_single[4,4] = 77.67
        ol_cij_single[5,5] = 78.36
        

        ol_cij = tex2elas.calc_cij(ol_cij_single, self.g_ol, 
                               self.volfrac_ol, scheme=scheme)
        return ol_cij
        
    def enstatite_cij(self, scheme='Voigt'):

        # Single crystal enstatite elasticity from DRex - should
        # check that this is up to date.
        en_cij_single = np.zeros((6,6))
        en_cij_single[0,0] = 236.9
        en_cij_single[1,1] = 180.5
        en_cij_single[2,2] = 230.4
        en_cij_single[0,1] = 79.6
        en_cij_single[1,0] = en_cij_single[0,1]
        en_cij_single[2,0] = 63.2
        en_cij_single[0,2] = en_cij_single[2,0]
        en_cij_single[1,2] = 56.8
        en_cij_single[2,1] = en_cij_single[1,2]
        en_cij_single[3,3] = 84.3
        en_cij_single[4,4] = 79.4
        en_cij_single[5,5] = 80.1
        

        en_cij = tex2elas.calc_cij(en_cij_single, self.g_en, 
                               self.volfrac_en, scheme=scheme)
        return en_cij

    def bulk_cij(self, scheme='Voigt'):

        if scheme=='Voigt':
            cij = self.fraction_olivine*self.olivine_cij(scheme='Voigt') + (1.0 - 
                 self.fraction_olivine)*self.enstatite_cij(scheme='Voigt')
        elif scheme == 'Reuss':
            cij = self.fraction_olivine*self.olivine_cij(scheme='Reuss') + (1.0 - 
                 self.fraction_olivine)*self.enstatite_cij(scheme='Reuss')
        elif scheme == 'Hill':
            cij = (self.bulk_cij(scheme='Voigt') + self.bulk_cij(scheme='Reuss'))/2.0
        else:
            raise ValueError('Scheme argument must be one of Voigt, Reuss or Hill')

        return cij

    def bulk_anisotropy(self, scheme='Voigt'):
        """Calculate and return the universal anisotropy index

        """
        Ua, Ua_sigma = CijUtil.uAniso(self.bulk_cij(scheme=scheme), 
                           np.zeros((6,6))) # We don't use the error (Ua_sigma) 
                                            # as we don't know the error on the
                                            # single crystal elasticity.

        return Ua

    def _unpack_drex_rdata(self, data, ngr):
        """Data is a 1D numpy array. This function pulls out the usefull info"""

        self.g_ol = np.copy(data[0:3*3*ngr].reshape((3,3,ngr)).T)
        self.g_en = np.copy(data[3*3*ngr:2*3*3*ngr].reshape((3,3,ngr)).T)
        self.volfrac_ol = np.copy(data[2*3*3*ngr:2*3*3*ngr+ngr])
        self.volfrac_en = np.copy(data[2*3*3*ngr+ngr:2*3*3*ngr+2*ngr])
        self.fraction_olivine = np.copy(data[3*3*ngr*2+ngr*4+10])

def process_drex_particles(particles):

    drex_particles = []
    for particle in particles:
        if particle['particle_class'] == 'drex':
            try:
                drex_particles.append(drex_particle(
                    particle['particle_position'].astype(np.float),
                    particle['particle_class'], particle['particle_rdata'].astype(np.float),
                    particle['particle_idata'].astype(np.int32)))
            except:
                print particle['particle_rdata']
                print "Skipped this particle"
                raise

    return drex_particles

def assert_is_rotmat(rotmat):
    """Checks that a numpy matrix is a rotation matrix

       Rotation matricies must be of shape (3,3), have a
       determinent of 1 and have the property of the 
       inverse being equal to the transpose. We check 
       all three and raise an AssertionError if this 
       is not the case."""

    assert rotmat.shape == (3,3)
    np.testing.assert_array_almost_equal(np.linalg.det(rotmat), 1.0)
    np.testing.assert_array_almost_equal(rotmat.transpose(), np.linalg.inv(rotmat))


def plot_particle_list(particle_list):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = []
    ys = []
    zs = []
    ua = []
    for particle in particle_list:
        xs.append(particle.position[0])
        ys.append(particle.position[1])
        zs.append(particle.position[2])
        ua.append(particle.bulk_anisotropy(scheme='Hill'))

    ax.scatter(xs, ys, zs, c=ua)
    plt.show()

if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser(description=
                'Read data from Theia particle data files, convert it to elasticity, and print to STDOUT')
    parser.add_argument('theiafile', help='Theia data file(s)', nargs='+')
    parser.add_argument('-s', '--scheme', help='Elasticity averaging scheme',
                choices=['Hill', 'Voigt', 'Reuss'], default='Hill')
    parser.add_argument('-p', '--plot', action='store_true', 
                help='Plot the universal anisotropy index distribution in 3D')
    #parser.add_argument('-o', '--outfile', help='Ouput plot to a file')
    args = parser.parse_args()

    # Read particles out of all files, concating output
    all_particles = []
    for this_file in args.theiafile:
        sys.stderr.write("Reading data from {0}\n".format(this_file))
        header_data, particles = read_texture_file(this_file)
        all_particles.extend(particles)

    sys.stderr.write("{0} particles read\n".format(len(all_particles)))
    drex_particles = process_drex_particles(all_particles)
    sys.stderr.write("{0} DRex particles processed\n".format(len(drex_particles)))

    if args.plot:
        plot_particle_list(drex_particles)


