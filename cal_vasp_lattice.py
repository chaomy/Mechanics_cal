#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_vasp_lattice.py
#
###################################################################
#
# Purpose :     multithreads to calculate lattice constant
#
# Creation Date :
# Last Modified : Thu Mar 30 23:56:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import numpy as np
import glob
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
from optparse import OptionParser
import os
import ase.io
import ase


try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import md_pot_data

except ImportError:
    print "error in importing"


class cal_lattice(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  gn_pbs.gn_pbs,
                  get_data.get_data,
                  gn_kpoints.gn_kpoints,
                  gn_incar.gn_incar):

    def __init__(self):
        self.figsize = (8, 6)
        self.npts = 20
        self.kpoints = [31, 31, 31]
        self._element = 'Nb'
        self.root = os.getcwd()
        self.pot = md_pot_data.dft_data.W_pbe

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        return

    def interpolate_lattice(self, filename):
        data = np.loadtxt(filename)

        lattice = data[:, 0]
        energy = data[:, 1]

        InterPoints = np.linspace(lattice[0], lattice[-1], 101)
        f = interp1d(lattice, energy)
        Ynew = f(InterPoints)

        i = np.argmin(Ynew)
        print "min energy ", np.min(energy)
        print "num", np.argmin(energy), "lattice", InterPoints[i]
        print (np.min(Ynew))

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        ax.plot(lattice, energy)
        plt.savefig("lattice.png")
        plt.show()
        return

    def loop_kpts(self):
        files = glob.glob('./DATA*')
        for file in files:
            cal_lattice(file)
        return

    def prepare_vasp_inputs(self, dirname):
        self.set_incar_type('dftunrelax')
        self.write_incar()

        self.set_diff_kpoints([31, 31, 31])
        self.set_intype('gamma')
        self.write_kpoints()

        self.set_pbs_type('va')
        self.set_wall_time(10)
        self.set_job_title(dirname)
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=None)
        os.system("cp ../../POTCAR .")
        return

    def gn_bcc(self):
        alat0 = self.pot['latbcc']
        delta = 0.001
        rng = [-15, 15]
        bcc_drv = gn_config.bcc(self.pot)

        for i in range(rng[0], rng[1]):
            alat = alat0 + i * delta

            if i >= 0:
                dirname = "dir-p-{:03d}".format(i)
            else:
                dirname = "dir-n-{:03d}".format(abs(i))

            self.mymkdir(dirname)
            os.chdir(dirname)

            bcc_drv.set_lattce_constant(alat)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')

            self.prepare_vasp_inputs(dirname)
            os.chdir(self.root)
        return

    def gn_fcc(self):
        alat0 = 3.90
        delta = 0.005
        #  rng = [-20, 20]
        rng = [20, 50]
        fcc_drv = gn_config.fcc(self.pot)

        for i in range(rng[0], rng[1]):
            alat = alat0 + i * delta
            if i >= 0:
                dirname = "dir-p-{:03d}".format(i)
            else:
                dirname = "dir-n-{:03d}".format(abs(i))

            self.mymkdir(dirname)
            os.chdir(dirname)

            fcc_drv.set_lattce_constant(alat)
            atoms = fcc_drv.set_fcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')

            self.prepare_vasp_inputs(dirname)
            os.chdir(self.root)
        return

    def collect_data(self, tag='hcp'):
        rng = [-15, 15]
        rng = [-8, 8]
        rng = [-20, 50]

        data = np.zeros([rng[1] - rng[0], 2])
        cnt = 0
        for i in range(rng[0], rng[1]):
            if i >= 0:
                dirname = "dir-p-{:03d}".format(i)
            else:
                dirname = "dir-n-{:03d}".format(abs(i))
            os.chdir(dirname)

            (energy, vol, atoms) = self.vasp_energy_stress_vol_quick()
            if (tag == 'fcc') or (tag == 'bcc'):
                data[cnt, 0] = atoms.get_cell()[0, 1]
                data[cnt, 1] = (energy)

            if tag == 'hcp':
                data[cnt, 0] = atoms.get_cell()[0, 0]
                data[cnt, 1] = (energy)

            cnt += 1
            os.chdir(self.root)
        np.savetxt('lat.dat', data)
        return

    def gn_hcp(self):
        alat0 = self.pot['ahcp']
        hcp_drv = gn_config.hcp(self.pot)
        for i in range(-10, 10):
            delta = 0.01
            alat = alat0 + i * delta
            if i >= 0:
                dirname = "dir-%03d" % (i)
            else:
                dirname = "dir-n%03d" % (i)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)

            os.chdir(dirname)
            hcp_drv.write_hcp_poscar(alat)

            self.set_incar_type("isif4")
            self.write_incar()

            self.set_diff_kpoints([23, 23, 23])
            self.set_intype('gamma')
            self.write_kpoints()

            self.set_pbs_type('va')
            self.set_wall_time(40)
            self.set_job_title(dirname)
            self.set_nnodes(1)
            self.set_ppn(12)
            self.set_main_job("mpirun vasp")
            self.write_pbs(od=None)

            os.system("cp ../../POTCAR .")
            os.chdir(self.root)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t",
                      "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      help="",
                      default="prp_r")

    (options, args) = parser.parse_args()
    drv = cal_lattice()
    if options.mtype.lower() == 'plt':
        drv.interpolate_lattice("lat.dat")

    elif options.mtype.lower() == 'loop':
        drv.loop_kpts()

    elif options.mtype.lower() == 'bcc':
        drv.gn_bcc()

    elif options.mtype.lower() == 'fcc':
        drv.gn_fcc()

    elif options.mtype.lower() == 'hcp':
        drv.gn_hcp()

    elif options.mtype.lower() == 'clcbcc':
        drv.collect_data('bcc')

    elif options.mtype.lower() == 'clchcp':
        drv.collect_data('hcp')

    elif options.mtype.lower() == 'clcfcc':
        drv.collect_data('fcc')
