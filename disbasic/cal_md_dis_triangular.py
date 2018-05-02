#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-05-01 22:36:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:38:45


class cal_md_triangular(object):

    def set_dipole_triangular_box(self, sizen=1):
        u = 1. / 3. * np.array([1., -2., 1.])
        v = 1. / 3. * np.array([2., -1., -1.])
        z = 1. / 3. * np.array([1., 1., 1.])

        # n1u = 3 * n - 1
        # n1v = 0
        # n2u = 0
        # n2v = 3 * n - 1
        # c1z = 1. / 3.
        # c2z = -c1z

        # c1 = n1u * u + n1v * v + c1z * z
        # c2 = n2u * u + n2v * v + c2z * z

        # print c1; print c2
        atoms = ase.lattice.cubic.BodyCenteredCubic(
            directions=[[1., -2., 1.],
                        [2., -1., -1.],
                        [1., 1., 1.]],
            latticeconstant=self.pot['latbcc'], size=(4, 4, 1),
            symbol=self.pot['element'])
        supercell = atoms.get_cell()
        addstrain = False
        if addstrain is True:
            supercell = strain * supercell
            atoms.set_cell(supercell)
            atoms.wrap(pbc=[1, 1, 1])
        ase.io.write('tri_perf_poscar.vasp', images=atoms, format='vasp')
        return atoms

    def bcc_screw_dipole_triangular_atoms(self, atoms=None, fname='qe.in'):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'],
            C12=self.pot['c12'],
            C44=self.pot['c44'])

        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        atoms = self.set_dipole_triangular_box()
        pos = atoms.get_positions()
        cell = atoms.get_cell()

        c1 = [0.51 * cell[0, 0], 1. / 3. * cell[1, 1]]
        c2 = [2 * cell[1, 0], 2. / 3. * cell[1, 1]]
        print(cell, c1, c2)

        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))
        ase.io.write('tri_dis_poscar.vasp', atoms, format='vasp')

    def cal_disp_dipo_lisa(self):  # never use it before
        e1 = 1. / 3. * np.array([-1., -1., 2.])
        e2 = 1. / 2. * np.array([1., -1., 0])
        e3 = 1. / 2. * np.array([1., 1., 1.])
        n, m = 21, 13

        v1 = n * e1 - 1. / (3. * m) * e3
        v2 = 0.5 * n * e1 + m * e2 + (0.5 - 1. / (6 * m)) * e3
        v3 = e3

        print(v1, v2, v3)
        v1 = np.round(v1, decimals=0)
        v2 = np.round(v2, decimals=0)
        v3 = np.round(v3, decimals=1)

        self.set_lattce_constant(self.pot['lattice'])
        self.set_element('Nb')
        print(v1, v2, v3)

        atoms = self.set_bcc_convention([v1, v2, v3], (n, 1, 1))
        ase.io.write("lmp_init.xyz", atoms, "xyz")
