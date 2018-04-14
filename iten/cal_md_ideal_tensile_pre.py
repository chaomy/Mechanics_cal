# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-07 13:15:10
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 14:02:00

import numpy as np
import ase
import ase.io


class iten_pre(object):

    def gn_primitive_lmps(self, strain=np.mat(np.identity(3)),
                          tag='lmp'):
        alat = self.alat
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])

        if tag == 'vasp':
            va_bas = bas * strain
            cell = alat * va_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR", images=atoms, format='vasp')

        if tag == 'qe':
            qe_bas = bas * strain
            cell = alat * qe_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.gn_qe_scf(atoms)

        if tag == 'lmp':
            lmp_bas = bas * strain
            lmp_bas = self.lmp_change_box(lmp_bas)
            cell = alat * lmp_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0., 0., 0.]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.write_lmp_config_data(atoms, 'init.txt')
