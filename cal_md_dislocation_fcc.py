#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-13 15:22:50


from ase import io
from ase.lattice import cubic
from numpy import array, sqrt
import atomman as am
import atomman.lammps as lmp


class md_dislocation_fcc(object):

    def fcc_screw_dipo(self):
        e1 = array([-2, -1., 1.])
        e2 = array([1., -1., 1.0])
        e3 = array([0.0, 1.0, 1.0])

        self.set_lattce_constant(4.05000)
        self.set_element('Al')

        atoms = cubic.FaceCenteredCubic(directions=[[-2, -1, 1],
                                                    [1, -1, 1],
                                                    [0, 1, 1]],
                                        size=(25, 10, 1),
                                        latticeconstant=4.050,
                                        symbol='Al',
                                        pbc=(1, 1, 1))

        movex = 0.0  # -np.sqrt(6.)/3. * self._alat
        io.write("lmp_perf.cfg",
                 atoms, "cfg")

        atoms = self.intro_dipole_screw_atoms(atoms)
        io.write("lmp_init.cfg",
                 atoms, "cfg")

        self.write_lmp_config_data(atoms)
        return

    def cal_cu3Au_dis(self):
        atoms = ase.io.read("./custom/perf.00000.dump",
                            format='lammps-dump')

        # delete some Cu atoms to create missing rows #
        unity = 3.63902984582228 * sqrt(2.) / 2.
        index_list = []
        for atom in atoms:
            if atom.position[1] < (8 * unity) and atom.symbol is 'H':
                index_list.append(atom.index)
        del atoms[index_list]

        # introduce edge dislocation #
        cell = atoms.get_cell()
        # zhigh = cell[2, 2]
        center = [0.5 * cell[0, 0], 10. * unity, 0.0]
        atoms = self.intro_single_edge_atoms(atoms,
                                             center, 1,
                                             [0, 1, 2])

        #  atoms = self.intro_edge_nuclei(atoms, center, -1,
        #  [0, 1, 2],
        #  [0.0*zhigh, zhigh])

        # cut extra half plane #
        #  atoms = self.cut_x_normal_atoms(atoms, self._alat,
        #  0, np.sqrt(2)/3.)

        io.write("ase.cfg", atoms, format='cfg')
        system, elements = am.convert.ase_Atoms.load(atoms)
        lmp.atom_data.dump(system, "lmp_init.txt")
        return
