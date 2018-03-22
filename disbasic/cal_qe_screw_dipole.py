#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-09 20:26:54


import numpy as np
import md_pot_data
import ase
import os
import gn_config
import get_data
import gn_pbs
import gn_qe_inputs
from . import cal_md_dis_dipole
from optparse import OptionParser


class qe_dislocation(get_data.get_data,
                     gn_qe_inputs.gn_qe_infile,
                     gn_pbs.gn_pbs, gn_config.bcc,
                     cal_md_dis_dipole.cal_dis_dipole):

    def __init__(self):
        self.pot = md_pot_data.qe_pot.vca_W75Re25
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        cal_md_dis_dipole.cal_dis_dipole.__init__(self, self.pot)
        return

    def gn_qe_screw_dipole_bcc(self):
        (dis_atoms, perf_atoms) = self.bcc_screw_dipole_configs_alongz()
        self.gn_infile_dipole_screw_atoms(dis_atoms)
        ase.io.write('dis_poscar', dis_atoms, format='vasp')
        ase.io.write('perf_poscar', perf_atoms, format='vasp')
        os.system("cp $POTDIR/{} . ".format(self.pot['file']))
        return

    def gn_infile_dipole_screw_atoms(self,
                                     atoms=None,
                                     fname='qe.in'):
        self.set_cal_type('relax')
        self.set_ecut('43')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-4')
        self.set_maxseconds(3600 * 80)
        with open(fname, 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons_tf(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (1, 2, 13))
            fid.close()
        return

    def cal_qe_dipo_screw(self,
                          in_tag="easy_easy",
                          input_s=0.0,
                          movex=0.0):
        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])
        n = 7
        m = 11
        atoms = self.set_bcc_convention([e1, e2, e3],
                                        (n, m, 1))  # z periodic 12
        atoms = self.cut_half_atoms(atoms)
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.5, 1.0, 0.5],
                         [0.0, 0.0, 1.0]])
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        atoms_perf = atoms.copy()
        atoms = self.intro_dipole_screw_atoms(atoms,
                                              self.pot['lattice'],
                                              movex,
                                              in_tag,
                                              input_s)
        self.write_poscar(atoms_perf)
        os.system("cp POSCAR new_perf.vasp")
        self.write_poscar(atoms)
        os.system("cp POSCAR POSCAR_init")
        self.gn_infile_dipole_screw_atoms(atoms)
        return (atoms_perf, atoms)

    def cal_qe_restart(self):
        atoms = self.qe_get_atom_pos_final()
        # atoms = self.qe_get_atom_pos_init()
        ase.io.write(filename='relax_poscar', images=atoms, format='vasp')
        # self.gn_infile_dipole_screw_atoms(atoms)
        return

    def cal_qe_restart_rescale_cell(self):
        atoms = self.qe_get_atom_pos()
        cell = atoms.get_cell()
        pot2 = md_pot_data.qe_pot.vca_W50Re50

        scaled_cell = cell / pot2['lattice']
        atoms.set_cell(scaled_cell * self.pot['lattice'])
        ase.io.write(filename='poscar_relax', images=atoms, format='vasp')
        # self.gn_infile_dipole_screw_atoms(atoms)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = qe_dislocation()

    dispatcher = {'dipole': drv.gn_qe_screw_dipole_bcc,
                  'restart': drv.cal_qe_restart,
                  're_cell': drv.cal_qe_restart_rescale_cell}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
