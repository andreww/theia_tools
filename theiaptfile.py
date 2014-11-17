#!/usr/bin/env python

import numpy as np

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
                        [line.rstrip('\r\n')[i:i+12] for i in xrange(0, len(line.rstrip('\r\n')), 12)]
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

    def _unpack_drex_rdata(self, data, ngr):
        """Data is a 1D numpy array. This function pulls out the usefull info"""

        self.g_ol = data[0:3*3*ngr].reshape((3,3,ngr)).T
        self.g_en = data[3*3*ngr:2*3*3*ngr].reshape((3,3,ngr)).T
        self.volfrac_ol = data[2*3*3*ngr:2*3*3*ngr+ngr]
        self.volfrac_en = data[2*3*3*ngr+ngr:2*3*3*ngr+2*ngr]

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
    for particle in particle_list:
        xs.append(particle.position[0])
        ys.append(particle.position[1])
        zs.append(particle.position[2])

    ax.scatter(xs, ys, zs)
    plt.show()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=
                'Read data from Theia particle data files.')
    parser.add_argument('theiafile', help='Theia data file')
    args = parser.parse_args()

    header_data, particles = read_texture_file(args.theiafile)
    print header_data
    print len(particles)
    drex_particles = process_drex_particles(particles)
    print len(drex_particles)
    plot_particle_list(drex_particles)

    print drex_particles[1].g_ol[0,:,:]
    print drex_particles[1].volfrac_ol[0]
    print np.sum(drex_particles[1].volfrac_ol[:])
    print drex_particles[1].g_en[0,:,:]
    print drex_particles[1].volfrac_en[0]
    print np.sum(drex_particles[1].volfrac_en[:])

    assert_is_rotmat(drex_particles[1].g_ol[0,:,:])
    assert_is_rotmat(drex_particles[1].g_en[10,:,:])
