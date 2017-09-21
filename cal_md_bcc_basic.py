#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-12 13:49:14


from optparse import OptionParser
from os import popen
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import md_pot_data
import plt_drv
import matplotlib.pyplot as plt

try:
    import atomman as am
    import atomman.lammps as lmp
    import gn_config
    import get_data
    import gn_lmp_infile
    import gn_pbs

except ImportError:
    print("error during import")


class cal_md_bcc_basic(gn_config.hcp,
                       gn_config.bcc,
                       gn_config.fcc,
                       get_data.get_data,
                       gn_pbs.gn_pbs,
                       gn_lmp_infile.gn_md_infile,
                       plt_drv.plt_drv):

    def __init__(self):
        gn_lmp_infile.gn_md_infile.__init__(self)
        self._element = 'Nb'
        self._pottype = 'eam/alloy'
        self._pot = '../Nb.eam.alloy.webarchive'
        self._lat = 3.308
        self.atoms = ase.lattice.cubic.BodyCenteredCubic(
            directions=[[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]],
            latticeconstant=self._lat,
            size=(1, 1, 1),
            symbol=self._element,
            pbc=(1, 1, 1))
        plt_drv.plt_drv.__init__(self)
        self.root = os.getcwd()
        return

    def loop_shear(self):
        for dim in range(20, 100, 10):
            atoms = ase.lattice.cubic.BodyCenteredCubic(
                directions=[[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
                latticeconstant=3.308,
                size=(10, dim, dim),
                symbol=self._element,
                pbc=(1, 1, 1))

            cell = atoms.get_cell()
            dirname = "dir-{}".format(dim)
            os.system("mkdir %s" % (dirname))

            for i in range(10):
                delta = i * 0.02
                strain = np.mat([[1.0, 0, 0.0],
                                 [delta, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]])
                new_cell = strain * cell

                atoms.set_cell(new_cell)

                pos = np.mat(atoms.get_positions())
                pos = pos * strain
                atoms.set_positions(pos)

                (system, elements) = am.convert.ase_Atoms.load(atoms)
                lmp.atom_data.dump(system, "lmp_init.txt")

                os.system("mpirun -n 4 lmp_mpi -i in.minimize")
                os.system("mv  bcc.dump  %s/dump_%03d" % (dirname, i))
            os.system("mv out.dat  out.dat.%03d" % (dim))
        return

    def cal_lattice(self, potname='pot.dat'):
        pot = md_pot_data.md_pot.Nb_adp
        popen("lmp_mpi -i in.bcc_adp")
        data = np.loadtxt("out.txt")
        pot['latbcc'], pot['ebcc'] = data[0], data[1]
        popen("lmp_mpi -i in.fcc_adp")
        data = np.loadtxt("out.txt")
        pot['latfcc'], pot['efcc'] = data[0], data[1]
        popen("lmp_mpi -i in.hcp_adp")
        data = np.loadtxt("out.txt")
        pot['ahcp'], pot['chcp'], pot['ehcp'] = data[0], data[1], data[2]

        print 'fcc', pot['efcc'] - pot['ebcc']
        print 'hcp', pot['ehcp'] - pot['ebcc']
        print 'latbcc', pot['latbcc']
        print 'latfcc', pot['latfcc']
        pot['lattice'] = pot['latbcc']
        self.pot = pot
        self.dump_data(potname, pot)
        return

    def loop_rcut_engy(self):
        npts = 16
        data = np.ndarray([npts, 2])
        for i in range(npts):
            rcut = 5.1500 + 0.01 * i
            dirname = 'dir-%5.4f' % (rcut)
            print dirname
            os.system("cp engy1000/{}/dummy.lamm*  .".format(dirname))
            potname = 'pot_%5.4f_lat' % (rcut)
            self.cal_lattice(potname)
            data[i, 0] = rcut
            data[i, 1] = self.pot['efcc'] - self.pot['ebcc']
        print data
        np.savetxt('rcut_ebcc2fcc.txt', data)
        return

    def plt_rcut_energy(self):
        data = np.loadtxt('rcut_ebcc2fcc.txt')
        plt.rc('xtick', labelsize=self.myfontsize - 2)
        plt.rc('ytick', labelsize=self.myfontsize - 2)
        self.set_keys()
        self.set_111plt()
        self.ax.plot(data[:, 0], data[:, 1],
                     **self.pltkwargs)
        self.ax.set_xlabel('cutoff radius (A)',
                           {'fontsize': self.myfontsize})
        self.ax.set_ylabel('fcc - bcc (meV / atom)',
                           {'fontsize': self.myfontsize})
        self.fig.savefig('rcut_ebcc2fcc.png', **self.figsave)
        return

    def cal_delta_energy(self):
        pot = md_pot_data.md_pot.Nb_adp
        print pot['efcc'] - pot['ebcc']
        print pot['ehcp'] - pot['ebcc']
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      help="",
                      default="prp_r")
    (options, args) = parser.parse_args()
    drv = cal_md_bcc_basic()
    if options.mtype.lower() == 'lattice':
        drv.cal_lattice()

    if options.mtype.lower() == 'shear':
        drv.loop_shear()

    if options.mtype.lower() == 'phasetrans':
        drv.cal_delta_energy()

    if options.mtype.lower() == 'looprcut':
        drv.loop_rcut_engy()

    if options.mtype.lower() == 'pltrcut':
        drv.plt_rcut_energy()
