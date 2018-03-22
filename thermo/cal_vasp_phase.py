#!/usr/bin/env python
# encoding: utf-8

import os
import glob
import numpy as np
from md_pot_data import dft_data
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

    def cal_bcc_to_fcc_struct(self):
        ilx, ily, ilz = 3, 3, 3

        pcub = []
        crit_dix = 1.5
        # [A]
        tworatio = 1.0

        sxx = 1.2
        syy = 1. / np.sqrt(sxx)

        strain = np.mat([[sxx, 0, 0],
                         [0, syy, 0],
                         [0, 0.0, syy]])

        pcub.append(np.array([0, 0, 0]) * self.lat)
        pcub.append(np.array([1, 0, 0]) * self.lat)
        pcub.append(np.array([0, 1, 0]) * self.lat)
        pcub.append(np.array([0, 0, 1]) * self.lat)
        pcub.append(np.array([1, 1, 0]) * self.lat)
        pcub.append(np.array([0, 1, 1]) * self.lat)
        pcub.append(np.array([1, 0, 1]) * self.lat)
        pcub.append(np.array([1, 1, 1]) * self.lat)

        while (True):
            atoms = ase.lattice.cubic.SimpleCubic(
                directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                size=(ilx, ily, ilz),
                symbol='Nb',
                pbc=(1, 1, 1),
                latticeconstant=self.lat)

            ########################################################
            # in Bcc you need to add 4 x 4 x 4 = 64
            # in Fcc you add 3 x 64;
            ########################################################

            for i in range(ilx):
                for j in range(ily):
                    for k in range(ilz):

                        if (np.random.rand() <= tworatio):
                            ### add two atoms ###
                            while (True):
                                add1 = np.random.rand(
                                    3) * self.lat
                                add2 = np.random.rand(
                                    3) * self.lat
                                a = []

                                for pos in pcub:
                                    a.append((
                                        np.linalg.norm(add1 - pos) > crit_dix))
                                    a.append((
                                        np.linalg.norm(add2 - pos) > crit_dix))

                                    a.append(
                                        np.linalg.norm(add1 - add2) > crit_dix)

                                if (all(np.array(a)) is True):
                                    break

                            add1[0] = (i) * self.lat + add1[0]
                            add1[1] = (j) * self.lat + add1[1]
                            add1[2] = (k) * self.lat + add1[2]

                            add2[0] = (i) * self.lat + add2[0]
                            add2[1] = (j) * self.lat + add2[1]
                            add2[2] = (k) * self.lat + add2[2]

                            atoms.append(
                                ase.Atom(
                                    symbol="Nb",
                                    position=(add1[0], add1[1], add1[2])))
                            atoms.append(
                                ase.Atom(
                                    symbol="Nb",
                                    position=(add2[0], add2[1], add2[2])))

                        else:
                            ### add one atom ###
                            while (True):
                                add1 = np.random.rand(
                                    3) * self.lat
                                a = []

                                for pos in pcub:
                                    a.append((
                                        np.linalg.norm(add1 - pos) > crit_dix))

                                    if (all(np.array(a)) is True):
                                        break

                            add1[0] = (i) * self.lat + add1[0]
                            add1[1] = (j) * self.lat + add1[1]
                            add1[2] = (k) * self.lat + add1[2]

                            atoms.append(
                                ase.Atom(symbol="Nb",
                                         position=(add1[0], add1[1], add1[2])))
                            a = []
                            b = []

            for atom1 in atoms[27:]:
                for atom2 in atoms[27:]:
                    diff = np.linalg.norm(atom1.position - atom2.position)

                    if (diff > 0.0000001):
                        b.append(diff)
                        a.append((diff > crit_dix))

            print(np.min(np.array(b)))
            if (all(np.array(a)) is True):
                break

        #### add strain ###
        #  cell = atoms.get_cell();
        #  cell = strain * cell;
        #  atoms.set_cell(cell, scale_atoms = True);
        atoms = self.volume_conserving_mono_strain_atoms(0.1, atoms)
        self.write_poscar(atoms)
        return

    def generate_rand_given_distance(self, r_min=2.0, r_max=2.3,
                                     dirname=None):
        ####### Spherical coordinates to Cartesian coordinates ########
        # theta = [0, 2pi);   phi = [0, pi)
        # x = r * cos theta * sin phi
        # y = r * sin theta * sin phi
        # z = r * cos phi
        ###############################################################
        # put a point in the center => put second => put third => ...

        #### initialize the atoms class and the supercell ####
        size = [4, 4, 4]
        lx, ly, lz = size[0] * self.lat, \
            size[1] * self.lat, \
            size[2] * self.lat

        supercell = np.mat([[lx, 0, 0], [0, ly, 0], [0, 0, lz]])

        pcub = []

        ##### all periodicities  6 share face      #####
        ##### 12 share length   and 8 shear corner #####
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

        for i in range((3 + 6 + 4)):
            pcub.append(-1 * pcub[i])

        atoms = ase.Atoms(cell=supercell, pbc=[1, 1, 1])
        print(atoms.get_cell())

        atom_num = 210
        m_pi = np.pi
        m_2pi = 2 * np.pi

        ### put the first atom in it ###
        atoms.append(
            ase.Atom(symbol="Nb", position=(0.5 * lx, 0.5 * ly, 0.5 * lz)))

        r_lim = 0.5 * r_min + 0.2

        if r_lim <= 2.0:
            r_lim = 2.0

        r_shl = r_max - r_min

        for i in range(atom_num - 1):
            cnt = 0
            while (True):
                theta = np.random.rand() * m_2pi
                phi = np.random.rand() * m_pi
                r = r_min + np.random.rand() * r_shl

                dx = r * cos(theta) * sin(phi)
                dy = r * sin(theta) * sin(phi)
                dz = r * cos(phi)

                pos_new = atoms[i - int(cnt / 30)].position + np.array(
                    [dx, dy, dz])
                cnt += 1
                if cnt >= 30 * i:
                    cnt = 0

                ###### add periodic  #######
                for k in range(3):
                    if (pos_new[k] > supercell[k, k]):
                        pos_new[k] -= supercell[k, k]
                    elif (pos_new[k] < 0):
                        pos_new[k] += supercell[k, k]

                ###### see whether larger than rmin #######
                a = []
                for atom in atoms:
                    a.append(np.linalg.norm(atom.position - pos_new) > r_lim)

                if (all(np.array(a)) is True):
                    #### check the distance with images  #####
                    b = []
                    for j in range(len(pcub)):
                        image = pcub[j]
                        a = []
                        for atom in atoms:
                            a.append(
                                np.linalg.norm(atom.position + image - pos_new)
                                > r_lim)
                            if (all(np.array(a)) is True):
                                b.append(True)
                        else:
                            b.append(False)

                    if ((len(b) == len(pcub)) and (all(np.array(b)) is True)):
                        break
                    ####  break the while loop ####

            atoms.append(ase.Atom(symbol="Nb", position=pos_new))
            print(len(atoms))
        self.write_poscar(atoms)
        self.prepare_dislocation_vasp_infiles(dirname)
        return
