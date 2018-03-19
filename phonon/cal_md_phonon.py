#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-07 13:09:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-17 08:56:13

import numpy as np
import ase
import ase.io
import os
import ase.lattice
import gn_config
from optparse import OptionParser
import md_pot_data
import atomman as am
import atomman.lammps as lmp
import get_data
import glob


class lmp_phonon(get_data.get_data):

    def __init__(self):
        get_data.get_data.__init__(self)
        self.pot = self.load_data("../pot.dat")

    def lmp_change_box(self, in_cell):
        unit_x = np.linalg.norm(in_cell[0, :])
        unit_y = np.linalg.norm(in_cell[1, :])
        unit_z = np.linalg.norm(in_cell[2, :])

        vect_x = in_cell[0, :] / unit_x
        vect_y = in_cell[1, :] / unit_y
        vect_z = in_cell[2, :] / unit_z

        ax = unit_x
        bx = np.dot(in_cell[1, :], vect_x.transpose())
        by = np.linalg.norm(np.cross(vect_x, in_cell[1, :]))

        cx = np.dot(in_cell[2, :],  vect_x.transpose())

        a_cross_b = np.cross(in_cell[0, :], in_cell[1, :])
        a_cross_b_unit = a_cross_b / np.linalg.norm(a_cross_b)
        Acvect_ycvect_x = np.cross(a_cross_b_unit, vect_x)

        cy = np.dot(in_cell[2, :], Acvect_ycvect_x.transpose())
        cz = np.abs(np.dot(in_cell[2, :], a_cross_b_unit.transpose()))

        out_cell = np.mat(np.zeros([3, 3], "float"))

        out_cell[0, 0], out_cell[0, 1], out_cell[0, 2] = ax, 0,  0
        out_cell[1, 0], out_cell[1, 1], out_cell[1, 2] = bx, by, 0
        out_cell[2, 0], out_cell[2, 1], out_cell[2, 2] = cx, cy, cz
        return out_cell

    def gn_primitive_lmps_with_strain(self):
        # 0.984257    0.983705
        strain = np.mat([[1.06, 0., 0.],
                         [0.0,  0.984257, 0.0],
                         [0.0, 0.0, 0.983705]])
        strain = np.mat([[1.10, 0., 0.],
                         [0.0,  0.974257, 0.0],
                         [0.0, 0.0, 0.973610]])
        self.gn_primitive_lmps(strain)

    def gn_primitive_lmps(self, strain=np.mat([[1., 0., 0.],
                                               [0., 1., 0.],
                                               [0., 0., 1.]])):
        alat = self.pot['lattice']
        cell = np.mat([[-0.5, 0.5, 0.5],
                       [0.5, -0.5, 0.5],
                       [0.5, 0.5, -0.5]])
        # [[-0.53       0.53       0.53     ]
        # [ 0.4921285 -0.4921285  0.4921285]
        # [ 0.4918525  0.4918525 -0.4918525]]
        cell = strain * cell
        print cell
        cell = alat * self.lmp_change_box(cell)

        atoms = ase.Atoms(
            'Nb', positions=[[0, 0, 0]], cell=cell, pbc=[1, 1, 1])
        ase.io.write("POSCAR", images=atoms, format='vasp')

        pos = np.array([[0, 0, 0]])

        file_name = 'lmp_init.txt'
        atom_num = 1
        with open(file_name, mode="w") as fout:
            fout.write("#lmp data config")
            fout.write("\n")
            fout.write("%d atoms\n" % (1))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (0, cell[0, 0]))
            fout.write("%f\t%f ylo yhi\n" % (0, cell[1, 1]))
            fout.write("%f\t%f zlo zhi\n" % (0, cell[2, 2]))
            fout.write("%f  %f  %f xy xz yz\n"
                       % (cell[1, 0],
                          cell[2, 0],
                          cell[2, 1]))

            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(atom_num):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                           % (i + 1, pos[i, 0], pos[i, 1], pos[i, 2]))
        fout.close()

#  Phonopy use the supercell base vector as basis
#  the basis vector we use is
#  [[-0.5, 0.5, 0.5],
#  [0.5, -0.5, 0.5],
#  [0.5, 0.5, -0.5]]
#  as we want to calculate
#  [[0.5, 0.5, 0.5]
#  [0.0, 0.0, 0.0],
#  [0.0, 0.5, 0.0],
#  [0.0, 0.5, 0.5]]

    def convert_kpath_upon_base_vector(self):
        unitcell = np.mat([[-0.5, 0.5, 0.5],
                           [0.5, -0.5, 0.5],
                           [0.5, 0.5, -0.5]])

        # add strain on it
        strain = np.mat([[1.06, 0., 0.],
                         [0.0,  0.984257, 0.0],
                         [0.0, 0.0, 0.983705]])

        # 0.974257    0.973610
        strain = np.mat([[1.10, 0., 0.],
                         [0.0,  0.974257, 0.0],
                         [0.0, 0.0, 0.973610]])
        strain = np.mat(np.identity(3))

        unitcell = strain * unitcell

        # for lammps we need to convert the unit cell
        #  unitcell = self.lmp_change_box(unitcell)
        # N[0.5 0.5 0.0]','\Gamma','H[1.0 0.0 0.0]','P[0.5 0.5 0.5]','\Gamma'
        # '\Gamma','H[1.0 0.0 0.0]','P[0.5 0.5 0.5]','\Gamma', N[0.5 0.5 0.0]',

        standardkpath = np.mat([[0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [0.5, 0.5, 0.5],
                                [0.0,  0., 0.],
                                [0.5, 0.5, 0.0]], 'float')

        transformedkpath = (unitcell * standardkpath.transpose()).transpose()
        stringlist = []
        for i in range(len(transformedkpath)):
            stringlist.append("[%3.3f, %3.3f, %3.3f]" %
                              (transformedkpath[i, 0],
                               transformedkpath[i, 1],
                               transformedkpath[i, 2]))
        print stringlist
        for i in range(len(stringlist) - 1):
            print "self.append_band(bands, %s, %s)" \
                % (stringlist[i],
                   stringlist[i + 1])

    def convert_pos_to_lmp_data_norm(self, filename="POSCAR-001"):
        ase_atoms = ase.io.read(filename, format='vasp')
        system, elements = am.convert.ase_Atoms.load(ase_atoms)
        lmp.atom_data.dump(system, "lmp_init.txt")

    def convert_pos_to_lmp_data(self, filename="POSCAR-001"):
        ase_atoms = ase.io.read(filename, format='vasp')
        cell = ase_atoms.get_cell()
        print np.linalg.norm(cell)
        print ase_atoms.get_positions()
        drv_gnconfig = gn_config.gnStructure(self.pot)
        drv_gnconfig.write_lmp_config_data(ase_atoms)

    def cal_phonon_auto(self):
        self.gn_primitive_lmps()
        os.system("phonopy --dim=\"5 5 5\" -d")
        os.system("lmp_mpi -i in.phonon")
        files = glob.glob("bcc.0.dump")
        os.system("phonopy -f {} --lammps".format(files[0]))
        os.system("MutiPhonon_oneAtom.py")


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = lmp_phonon()
    dispatcher = {'cnv': drv.convert_pos_to_lmp_data,
                  'prim': drv.gn_primitive_lmps,
                  'strain': drv.gn_primitive_lmps_with_strain,
                  'kpt': drv.convert_kpath_upon_base_vector,
                  'auto': drv.cal_phonon_auto}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
