#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-23 22:22:14


import os
import glob
import numpy as np
from md_pot_data import va_pot
import ase
from optparse import OptionParser

try:
    import gn_config
    import get_data
    import gn_incar
    import gn_kpoints
    import gn_pbs
    import Intro_vasp

except ImportError:
    print("error during import")


class vasp_inters(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  get_data.get_data,
                  gn_incar.gn_incar,
                  gn_pbs.gn_pbs,
                  gn_kpoints.gn_kpoints,
                  Intro_vasp.vasp_change_box):

    def __init__(self, structure=None):
        self.pot = va_pot.W_pbe
        self.pot['lattice'] = self.pot['latbcc']
        gn_incar.gn_incar.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)

        e1 = np.array([1., 0., 0.]) * self.pot['lattice']
        e2 = np.array([0., 1., 0.]) * self.pot['lattice']
        e3 = np.array([0., 0., 1.]) * self.pot['lattice']
        #  self.set_hcp_lat(self.pot['lattice'],
        #  np.sqrt(8) / 3. * self.pot['lattice'])
        self.root = os.getcwd()

        #############################################################
        # self.atoms is body centered cubic
        # self.atoms_new is simple cubic
        #############################################################
        self.atoms = self.set_bcc_convention([[e1[0], e1[1], e1[2]],
                                              [e2[0], e2[1], e2[2]],
                                              [e3[0], e3[1], e3[2]]],
                                             (1, 1, 1))

        atoms_new = self.atoms.copy()
        for atom in self.atoms:
            if atom.position[1] > 1.0:
                print atom.index
                del atoms_new[atom.index]
        self.atoms_simple = atoms_new
        return

    def prepare_dislocation_vasp_infiles(self, dirname='inters'):
        self.set_incar_type("dftrelax")
        self.write_incar()

        self.set_diff_kpoints([4, 4, 4])
        self.write_kpoints()

        self.set_nnodes(4)
        self.set_pbs_type('va')
        self.set_wall_time(140)
        self.set_job_title(dirname)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=False)

        os.system("cp ../POTCAR  .")
        return

    def cal_dumbbell_100(self, atoms):
        alat = self.pot['lattice']
        dirname = "dir-dumb100"
        self.mymkdir(dirname)
        os.chdir(dirname)

        #  x = (np.sqrt(10) - 1)/3.
        #  b = 0.5 * (1 - x)
        b = 0.25

        pos1 = np.array([0.5, 0.5, b]) * alat
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        pos2 = np.array([0.5, 0.5, 1 - b]) * alat
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos2))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_dumbel_100.vasp")

        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def cal_dumbbell_110(self, atoms):
        dirname = "dir-dumb110"
        self.mymkdir(dirname)
        os.chdir(dirname)

        b1 = (4 * np.sqrt(2) - np.sqrt(11)) / (6 * np.sqrt(2))
        b2 = 1 - b1

        # add atoms #
        pos1 = np.array([b1, b1, 0.5]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        pos2 = np.array([b2, b2, 0.5]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos2))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_dumbel_110.vasp")

        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def cal_dumbbell_111(self, atoms):
        dirname = "dir-dumb111"
        self.mymkdir(dirname)

        os.chdir(dirname)

        b1 = 1. / 3.
        b2 = 2. / 3.

        # add atoms #
        pos1 = np.array([b1, b1, b1]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        pos2 = np.array([b2, b2, b2]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos2))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_dumbel_111.vasp")
        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def cal_crowdion(self, atoms):
        dirname = "dir-crowdion"
        self.mymkdir(dirname)
        os.chdir(dirname)

        b1 = 0.25
        # add atoms #
        pos1 = np.array([b1, b1, b1]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_crowdion.vasp")

        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def cal_octahedral(self, atoms):
        dirname = "dir-octahedral"
        self.mymkdir(dirname)
        os.chdir(dirname)

        b = 0.5
        pos1 = np.array([b, b, 0]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_octahedral.vasp")

        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def cal_tetrahedral(self, atoms):
        dirname = "dir-tetrahedral"
        self.mymkdir(dirname)
        os.chdir(dirname)

        pos1 = np.array([0.25, 0.5, 0]) * self.pot['lattice']
        atoms.append(ase.atom.Atom(self._element,
                                   position=pos1))

        self.write_poscar(atoms)
        os.system("cp POSCAR ../POSCAR_tetrahedral.vasp")

        self.prepare_dislocation_vasp_infiles(dirname)
        os.chdir(self.root)
        return

    def surf(self, tag):
        if tag == "read":
            atoms = ase.io.read("POSCAR", format='vasp')

            #  atoms = atoms.repeat((3, 3, 1))  # surface 100
            #  atoms = atoms.repeat((2, 2, 1))  # surface 110
            atoms = atoms.repeat((1, 2, 1))  # surface 111

        # add perturbation #
        cnt = 0
        while True:
            atoms_copy = self.add_perturbation(atoms, 0.3)  # 0.5;  0.4
            if self.check_min(atoms_copy):
                atoms = atoms_copy
                break
            elif cnt > 20:
                print "no write"
                break
            else:
                cnt += 1

        ase.io.write(filename="POSCAR.vasp", images=atoms, format='vasp')
        return

    def monos(self, tag='bcc'):
        self.pot['latbcc'] = 0.85**(1. / 3.) * self.pot['latbcc']
        #  strain = [1.26, 1.0196, 0.8153]
        strain = [1., 1., 1.]
        # e1 = np.array([1., 0., 0.]) * lat * strain[0]
        # e2 = np.array([0., 1., 0.]) * lat * strain[1]
        # e3 = np.array([0., 0., 1.]) * lat * strain[2]
        if tag in ["bcc"]:
            atoms = self.set_bcc_convention([[1, 0, 0],
                                             [0, 1, 0],
                                             [0, 0, 1]], (4, 4, 4))
        elif tag in ["fcc"]:
            self.set_lattce_constant(self.pot['latfcc'])
            atoms = self.set_fcc_convention([e1, e2, e3], (3, 3, 3))
        elif tag in ["hcp"]:
            self.set_hcp_lattice_constant(self.pot['ahcp'],
                                          self.pot['chcp'])
            atoms = self.set_hcp_convention((4, 4, 4))
        elif tag == "read":
            atoms = ase.io.read("POSCAR", format='vasp')
            style = 'tp'
            if style == 'op':
                atoms = atoms.repeat((4, 3, 3))  # OP
            elif style == 'tp':
                atoms = atoms.repeat((4, 4, 4))  # TP
            elif style == 'shear':
                atoms = atoms.repeat((5, 5, 5))  # shear

        elif tag == "prim":
            atoms = self.set_bcc_primitive((5, 5, 5))
        elif tag == "vac":
            atoms = self.set_bcc_convention([e1, e2, e3], (5, 5, 5))
            atoms.pop(150)

        #  add strain
        #  atoms = self.volume_conserving_ortho_strain_atoms(0.03, atoms)
        #  atoms = self.volume_conserving_mono_strain_atoms(0.03, atoms)
        ###################################################################
        # add perturbation
        ###################################################################
        cnt = 0
        print atoms.get_positions()

        add_perturbation = True
        if add_perturbation is True:
            while True:
                atoms_copy = self.add_perturbation(atoms, 0.1)  # 0.5;  0.4
                if self.check_min(atoms_copy):
                    atoms = atoms_copy
                    break
                elif cnt > 20:
                    print "no write"
                    break
                else:
                    cnt += 1
        ase.io.write(filename="POSCAR.vasp",
                     images=atoms,
                     format='vasp')
        return

    def collect_energy(self):
        dir_list = glob.glob("dir-*")
        for dirname in dir_list:
            os.chdir(dirname)

            (energy, vol, atomnum) = self.vasp_energy_stress_vol_quick()
            print dirname
            print(energy / atomnum)
            print
            os.chdir(self.root)
        return

    def check_min(self, in_atoms=None):
        if in_atoms is None:
            atoms = ase.io.read("POSCAR_cnt", format="vasp")
        else:
            atoms = in_atoms
        distances = atoms.get_all_distances()
        minv = 10
        for i in range(len(distances)):
            for j in range(len(distances[i])):
                if (distances[i][j] < minv) and (distances[i][j] != 0):
                    minv = distances[i][j]
        if minv >= 2.0:
            return True
        else:
            return False

    def check_min_bad(self, in_atoms=None):
        if in_atoms is None:
            atoms = ase.io.read("POSCAR_cnt", format="vasp")
        else:
            atoms = in_atoms

        cell = atoms.get_cell()
        lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]

        pcub = []
        pcub.append(np.array([lx, 0, 0]))
        pcub.append(np.array([0, ly, 0]))
        pcub.append(np.array([0, 0, lz]))

        pcub.append(np.array([lx, ly, 0]))
        pcub.append(np.array([lx, 0, lz]))
        pcub.append(np.array([0, ly, lz]))

        pcub.append(np.array([lx, -ly, 0]))
        pcub.append(np.array([lx, 0, -lz]))
        pcub.append(np.array([0, ly, -lz]))

        pcub.append(np.array([lx, ly, lz]))
        pcub.append(np.array([-lx, ly, lz]))
        pcub.append(np.array([lx, -ly, lz]))
        pcub.append(np.array([lx, ly, -lz]))

        print("start check min")
        for atom1 in atoms:
            for atom2 in atoms:
                for image in pcub:
                    dis = np.linalg.norm(atom1.position - atom2.position +
                                         image)
                    if (dis > 0 and dis < 2.0):
                        print dis
                        return False
                    dis = np.linalg.norm(atom1.position - atom2.position -
                                         image)
                    if (dis > 0 and dis < 2.0):
                        print dis
                        return False
        return True

    def prep_interstitials(self):
        # reported size 5 x 5 x 5
        alat = self.pot['lattice']
        sizen = 5
        atoms = self.atoms.copy().repeat((sizen, sizen, sizen))
        pos = np.array([alat,  alat,  alat]) * 0.5
        for atom in atoms:
            if ((pos - atom.position) == ([0, 0, 0])).all():
                index = atom.index
                print "del {}  atom".format(index)
        del atoms[index]

        ase.io.write("POSCAR_vac", images=atoms, format='vasp')

        self.cal_dumbbell_100(atoms.copy())
        self.cal_dumbbell_110(atoms.copy())
        self.cal_dumbbell_111(atoms.copy())

        atomsper = self.atoms.copy().repeat((sizen, sizen, sizen))

        self.cal_crowdion(atomsper.copy())
        self.cal_octahedral(atomsper.copy())
        self.cal_tetrahedral(atomsper.copy())
        return


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

if __name__ == '__main__':
    drv = vasp_inters()

    if options.mtype.lower() == 'prep':
        drv.prep_interstitials()

    if options.mtype.lower() == "point_get":
        drv.collect_energy()

    if options.mtype.lower() == "collect":
        drv.loop_collect_r()

    if options.mtype.lower() == "checkmin":
        drv.check_min()

    if options.mtype.lower() == "monos":
        drv.monos()

    if options.mtype.lower() == "surf":
        drv.surf('read')
