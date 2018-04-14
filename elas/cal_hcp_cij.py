#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-09 18:51:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-13 20:44:01

from optparse import OptionParser
import os
import ase.io
import ase
import numpy as np
import gn_pbs
import gn_config
import md_pot_data
import ase.lattice.orthorhombic as otho
import ase.lattice.hexagonal as Hexagonal


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 1. / 3., 0.5],
                     [0.5, 5. / 6., 0.5]]

# class othoHCPFractory(otho.SimpleOrthorhombicFactory):
#     bravais_basis = [[0.0, 0.0, 0.0],
#                      [0.5, 0.5, 0.0],
#                      [1. / 3., 0.0, 0.5],
#                      [5. / 6., 0.5, 0.5]]

othoHCP = othoHCPFractory()
evA3toGpa = 160.21766208


class cal_cij_hcp(gn_config.gnStructure, gn_pbs.gn_pbs):

    def __init__(self):
        self.pot = md_pot_data.md_pot.mg_Poco
        gn_pbs.gn_pbs.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)

    def gn_hcp(self):
        aa, cc = self.pot['ahcp'], self.pot['chcp']
        atoms = Hexagonal.HexagonalClosedPacked(
            latticeconstant={'a': aa, 'c': cc},
            size=(1, 1, 1), symbol=self.pot["element"], pbc=(1, 1, 1))
        return atoms

    def gn_otho_hcp(self):
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'],
                                         np.sqrt(3) * self.pot['ahcp'],
                                         self.pot['chcp']),
                        size=(1, 1, 1), symbol=self.pot['element'])
        # atoms = othoHCP(latticeconstant=(np.sqrt(3) * self.pot['ahcp'],
        #                                  self.pot['ahcp'],
        #                                  self.pot['chcp']),
        #                 size=(1, 1, 1), symbol=self.pot['element'])
        return atoms

    def add_strain(self, kk, d, atoms):
        strains = {'B': np.mat([[1 + d, 0, 0],
                                [0, 1 + d, 0],
                                [0, 0, 1 + d]]),  # (2C11+2C12+C33+4C13)/2.
                   'Ca': np.mat([[1. + d, 0, 0],
                                 [0, 1 - d, 0.],
                                 [0, 0, 1.]]),   # C11 - C12
                   'Cb': np.mat([[1.0, 0.0, 0.0],
                                 [0.0, 1.0, 0.0],
                                 [d, 0.0, 1.0]]),  # C44 / 2.
                   'Cc': np.mat([[1 + d, 0, 0],
                                 [0, 1 + d, 0],
                                 [0, 0, 1 - 2 * d]]),  # (C11+C12+2C33-4C13)
                   'D': np.mat([[1, 0, 0],
                                [0, 1., 0.],
                                [0, 0, 1 + d]])}  # C33 / 2.
        strain = strains[kk]
        org_cell = atoms.get_cell()
        pos = np.mat(atoms.get_positions())
        # pos = pos * strain
        atoms.set_cell(strain * org_cell, scale_atoms=True)
        # atoms.set_positions(pos)

    def set_pbs(self, mdir):
        self.set_pbs_type('va')
        self.set_wall_time(8)
        self.set_job_title(mdir)
        self.set_nnodes(1)
        self.set_ppn(8)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=True)

    def prepare_vasp_inputs(self, mdir):
        self.set_pbs(mdir)
        os.system("mv POSCAR {}".format(mdir))
        os.system("mv va.pbs {}".format(mdir))
        os.system("cp ../KPOINTS {}".format(mdir))
        os.system("cp ../INCAR {}".format(mdir))
        os.system("cp ../POTCAR {}".format(mdir))

    def va_prep(self):
        for kk in ['B', 'Ca', 'Cb', 'Cc', 'D']:
            data = np.ndarray([10, 2])
            mdir = 'DIR_{}'.format(kk)
            self.mymkdir(mdir)
            os.chdir(mdir)
            for i in range(10):
                d = i * 0.0005 - 5 * 0.0005
                atoms = self.gn_hcp()
                self.add_strain(kk, d, atoms)
                ase.io.write(filename="POSCAR", images=atoms, format='vasp')
                tmp = "dir_{:03d}".format(i)
                self.mymkdir(tmp)
                self.prepare_vasp_inputs(tmp)
            os.chdir(os.pardir)

    def md_cals(self):
        for kk in ['B', 'Ca', 'Cb', 'Cc', 'D']:
            data = np.ndarray([10, 2])
            mdir = 'dir_{}'.format(kk)
            self.mymkdir(mdir)
            os.system("cp in.init {}".format(mdir))
            os.chdir(mdir)
            for i in range(10):
                d = i * 0.0005 - 5 * 0.0005
                # atoms = self.gn_otho_hcp()
                atoms = self.gn_otho_hcp()
                self.add_strain(kk, d, atoms)
                self.write_lmp_config_data(atoms)
                os.system('lmp_mpi -i in.init')
                data[i, :] = d, np.loadtxt("out.txt")
            np.savetxt("data.txt", data, fmt="%.10f")
            os.chdir(os.pardir)

    def fit(self):
        yy = np.zeros(5)
        vol = np.linalg.det(self.gn_otho_hcp().get_cell())
        for kk, i in zip(['B', 'Ca', 'Cb', 'Cc', 'D'], range(5)):
            data = np.loadtxt("dir_{}/data.txt".format(kk))
            data[:, 1] /= vol
            res = np.polyfit(data[:, 0], data[:, 1], deg=2)
            print(res)
            yy[i] = res[0]
        yy[0] *= 2.0   # (2C11 + 2C12 + C33 + 4C13)/2.
        yy[1] *= 1.0   # C11 - C12
        yy[2] *= 2.0   # C44 / 2.
        yy[3] *= 1.0   # (C11 + C12 + 2C33 - 4C13)
        yy[4] *= 2.0   # C33 / 2.
        print("C11, C12, C33, C13, C44")
        A = np.mat([[2., 2., 1., 4., 0.],
                    [1., -1., 0., 0., 0.],
                    [0., 0., 0., 0., 1.],
                    [1., 1., 2., -4., 0.],
                    [0., 0., 1., 0., 0]])
        cij = np.linalg.inv(A) * np.mat(yy).transpose()
        print(cij.transpose() * evA3toGpa)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_cij_hcp()
    dispatcher = {'cal': drv.md_cals,
                  'fit': drv.fit,
                  'va': drv.va_prep}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
