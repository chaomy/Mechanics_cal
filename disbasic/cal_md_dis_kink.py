#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-13 23:59:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:55:21


import atomman as am
import atomman.lammps as lmp
from . import cal_md_dis_schmid
import get_data
import gn_config
import Intro_vasp
import gn_pbs
import gn_lmp_infile
import os
import numpy as np
from optparse import OptionParser
import ase.io
import ase
import md_pot_data
import glob


class bcc_kink(gn_config.bcc,
               get_data.get_data,
               gn_pbs.gn_pbs,
               Intro_vasp.vasp_change_box,
               gn_lmp_infile.gn_md_infile):

    def __init__(self, structure='bcc'):
        gn_pbs.gn_pbs.__init__(self)
        self.pot = md_pot_data.md_pot.Nb_eam
        self._alat = self.pot['lattice']  # Nb
        Intro_vasp.vasp_change_box.__init__(self, self._alat)
        gn_lmp_infile.gn_md_infile.__init__(self)
        self.shmid_drv = cal_md_dis_schmid.cal_bcc_schmid(self.pot)

    def intro_kink_pair(self):
        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        atoms = self.set_bcc_convention([e1, e2, e3], (30, 30, 60))
        xc1 = (0.0 + (-2.56656)) / 2. + 45 * \
            np.sqrt(6.) / 3. * self.pot['lattice']
        yc1 = (0.0 + (2.22271)) / 2. + 15 * np.sqrt(2.) * self.pot['lattice']
        H = np.sqrt(2. / 3.0) * self.pot['lattice']

        h = 0.0 * H
        atoms = self.intro_kink_screw_dislocations(
            atoms, (xc1, yc1), (xc1 + H, yc1), h, 1. / 4.)

        ase.io.write("lmp_init.cfg", atoms, "cfg")
        fname = "init.data"
        self.write_lmp_config_data(atoms, fname)
        self.gn_md_minimize_cfg("init.data", "./w_eam4.fs", "W")

        ############################################################
        # cluster method  (large supercell) # make a plate # run some MD
        ############################################################
    def construct_init_final(self):
        alat = self._alat
        move1 = [-alat * np.sqrt(6.) / 3., 0, 0]
        self.shmid_drv.make_screw_plate(size=[60, 84, 20],
                                        rad=[152, 160],
                                        move=move1,
                                        tag='[211]',
                                        filename="init.txt")
        move2 = [0, 0, 0]
        self.shmid_drv.make_screw_plate(size=[60, 84, 20],
                                        rad=[152, 160],
                                        move=move2,
                                        tag='[211]',
                                        filename="final.txt")

        # dipole method
    def cal_kink_pair_neb_pre(self):
        dirname = "kink_0.004"
        self.mymkdir(dirname)
        os.chdir(dirname)

        e1 = 1. / 3. * np.array([1.,  1., -2.])
        e2 = 1. / 2. * np.array([-1., 1.,  0])
        e3 = np.array([0.5,  0.5,  0.5])

        Lx = 30  # Default 30
        Ly = 15  # Default 30
        Lz = 100

        xc1 = (0.0 + (-2.56656)) / 2. + 1.5 * \
            Lx * np.sqrt(6.) / 3. * self._alat
        yc1 = (0.0 + (2.22271)) / 2. + 15 * np.sqrt(2.) * self._alat
        H = np.sqrt(2. / 3.0) * self._alat

        npt1 = 8
        npt2 = 16
        delta = 1. / npt1

        neb_tag = "each"
        if (neb_tag == "each"):
            for i in range(npt2):
                #  first 14 configs  #
                atoms = self.set_bcc_convention([e1, e2, e3], (Lx, Ly, Lz))

                if (i == 0):
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               0)
                elif (i <= npt1):
                    h = 0.5 * H * i * delta
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               h, 1. / 3.)
                elif (i < (npt2 - 1)):
                    h = 0.5 * H
                    tilpnt = 1. / (3. + i - npt1)
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               h, tilpnt)
                else:
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1 + H, yc1),
                                                               (xc1 + 2 * H, yc1),
                                                               0)

                # Be careful which file lammps reads #
                fname = "init_%d.data" % (i + 1)
                coordn = "coord_%d" % (i + 1)
                self.write_lmp_config_data(atoms, file_name=fname)
                self.write_lmp_coords(atoms, file_name=coordn)

        elif (neb_tag == "final"):
            #  generate initial #
            atoms = self.set_bcc_convention([e1, e2, e3], (Lx, Ly, Lz))
            atoms = self.intro_kink_screw_dislocations(
                atoms, (xc1, yc1), (xc1 + H, yc1), 0)
            fname = "init.data"
            self.write_lmp_config_data(atoms, file_name=fname)

            # generate final #  atoms = self.set_bcc_convention([e1, e2, e3],
            # (30, 30, 100));
            atoms = self.intro_kink_screw_dislocations(
                atoms, (xc1 + H, yc1), (xc1 + 2 * H, yc1), 0)
            fname = "final.data"
            coordn = "final.screw"
            self.write_lmp_config_data(atoms, file_name=fname)
            self.write_lmp_coords(atoms, file_name=coordn)
        else:
            h = 0.5 * H
            atoms = self.set_bcc_convention([e1, e2, e3], (Lx, Ly, Lz))
            atoms = self.intro_kink_screw_dislocations(
                atoms, (xc1, yc1), (xc1 + H, yc1), h, 1. / 4.)
            fname = "init.data"
            self.write_lmp_config_data(atoms, file_name=fname)

        ############################################################
        # prepare for neb calculation of kink pair (disocation dipole)
        ############################################################
    def cal_kink_pair_neb_dipole(self, tag='final'):
        torient = 'y'
        if torient == 'y':
            e1 = 1. / 3. * np.array([1., 1., -2.])
            e2 = np.array([0.5, 0.5, 0.5])
            e3 = 1. / 2. * np.array([1, -1, 0])

        times = 4
        n = 7 * times
        m = 11 * times
        t_thick = 50

        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[e1, e2, e3],
                                                    latticeconstant=self._alat,
                                                    size=(n, t_thick,  m),
                                                    symbol=self._element,
                                                    pbc=(1, 1, 1))

        ###  cut the supercell (only need half) ###
        atoms = self.cut_half_atoms_new(atoms, "cutz")

        ###  add strain to cancell the strain field ###
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.5, 0.0, 1.0]])

        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])

        H = np.sqrt(2. / 3.0) * self._alat
        h = 0.5 * H
        tilpnt = 1. / 5.

        lxid = 0  # glide plane
        lyid = 2  # y
        lzid = 1  # dislocation line

        add_disp = [0.25 * (supercell[lxid, lxid] + supercell[lyid, lxid]),
                    0.50 * (supercell[lxid, lyid] + supercell[lyid, lyid])]

        center_opt = 'cal'
        if center_opt == 'box':
            xc1 = 0 + add_disp[0]
            yc1 = 0 + add_disp[1]

            xc2 = 0.5 * supercell[lxid, lxid] + add_disp[0]
            yc2 = 0.0 * supercell[lyid, lyid] + add_disp[1]

        elif center_opt == 'cal':
            unit_a = H
            unit_b = np.sqrt(2.) / 2. * self._alat
            print(unit_a, unit_b)
            yc1 = yc2 = (m / 2 + 0.5) * unit_b
            xc1 = (10 * times + 0.5) * unit_a
            xc2 = (10 * times + 0.5 + 1.5 * n) * unit_a
            print(xc1, yc1, xc2,  yc2)

        if tag == 'initial':
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc1, yc1),
                                                       (xc1 + H, yc1),
                                                       h, tilpnt, 1)

            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc2, yc2),
                                                       (xc2 + H, yc2),
                                                       h, tilpnt, -1)
            self.write_lmp_config_data(atoms, 'init.data')
            os.system("mv init.data  output")

        elif tag == 'final':
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc1 + H, yc1),
                                                       (xc1 + 2 * H, yc1),
                                                       h, tilpnt, 1)
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc2 + H, yc2),
                                                       (xc2 + 2 * H, yc2),
                                                       h, tilpnt, -1)
            #  self.write_lmp_config_data(atoms, 'final.txt')
            self.write_lmp_coords(atoms, file_name='final.coord')
            os.system("mv  final.coord  output")

        elif tag == "each_interp":
            atoms_copy = atoms.copy()

            ### before move ###
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc1, yc1),
                                                       (xc1, yc1),
                                                       0, tilpnt, 1)

            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc2, yc2),
                                                       (xc2, yc2),
                                                       0, tilpnt, -1)
            atoms_init = atoms.copy()
            atoms = atoms_copy.copy()

            ### after move ###
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc1 + H, yc1),
                                                       (xc1 + H, yc1),
                                                       0, tilpnt, 1)

            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc2 + H, yc2),
                                                       (xc2 + H, yc2),
                                                       0, tilpnt, -1)

            atoms_final = atoms.copy()
            atoms = atoms_copy.copy()

            ### during move ###
            h = 0.5 * H
            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc1, yc1),
                                                       (xc1 + H, yc1),
                                                       h, tilpnt, 1)

            atoms = self.intro_kink_screw_dislocations(atoms,
                                                       (xc2, yc2),
                                                       (xc2 + H, yc2),
                                                       h, tilpnt, -1)

            atoms_mid = atoms.copy()
            atoms = atoms_copy.copy()

            n1 = 8
            n2 = 8

            ##### init to transition state ####
            init_pos = atoms_init.get_positions()
            delta = atoms_mid.get_positions() - atoms_init.get_positions()
            delta = delta / n1

            for i in range(n1):
                pos_interp = init_pos + i * delta

                atoms.set_positions(pos_interp)

                fname = "init_%d.data" % (i + 1)
                coordn = "coord_%d" % (i + 1)
                self.write_lmp_config_data(atoms,
                                           file_name=fname)
                self.write_lmp_coords(atoms,
                                      file_name=coordn)

            ##### transition to final ####
            init_pos = atoms_mid.get_positions()
            delta = atoms_final.get_positions() - atoms_mid.get_positions()
            delta = delta / (n2 - 1)

            for i in range(n2):
                pos_interp = init_pos + i * delta

                atoms.set_positions(pos_interp)

                fname = "init_%d.data" % (n1 + i + 1)
                coordn = "coord_%d" % (n1 + i + 1)
                self.write_lmp_config_data(atoms, file_name=fname)
                self.write_lmp_coords(atoms, file_name=coordn)

        elif tag == "each_change_h":
            # this method is to change the h and tilpnt to generate the
            # configuration of each time of calculation
            dirname = "kinkPPPS6"
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            os.chdir(dirname)

            npt1 = 8
            npt2 = 24
            delta = 1. / npt1
            atoms_copy = atoms.copy()
            for i in range(npt2):
                atoms = atoms_copy.copy()
                if (i == 0):
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               0, tilpnt, 1)

                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc2, yc2),
                                                               (xc2 + H, yc2),
                                                               0, tilpnt, -1)

                elif (i <= npt1):
                    h = 0.5 * H * i * delta
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               h, tilpnt, 1)

                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc2, yc2),
                                                               (xc2 + H, yc2),
                                                               h, tilpnt, -1)

                elif (i < (npt2 - 1)):
                    h = 0.5 * H
                    tilpnt = 1. / (3. + i - npt1)
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1, yc1),
                                                               (xc1 + H, yc1),
                                                               h, tilpnt, 1)

                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc2, yc2),
                                                               (xc2 + H, yc2),
                                                               h, tilpnt, -1)

                else:
                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc1 + H, yc1),
                                                               (xc1 + H, yc1),
                                                               0, tilpnt, 1)

                    atoms = self.intro_kink_screw_dislocations(atoms,
                                                               (xc2 + H, yc2),
                                                               (xc2 + H, yc2),
                                                               0, tilpnt, -1)

                ###### Be careful which file lammps reads ######
                fname = "init_%d.data" % (i + 1)
                coordn = "coord_%d" % (i + 1)
                self.write_lmp_config_data(atoms, file_name=fname)
                self.write_lmp_coords(atoms, file_name=coordn)

        #  elif tag == 'perf':
            #  self.write_lmp_config_data(atoms, 'lmp_init.txt')

    def cal_kink_energy(self):
        atoms = ase.io.read('./relax.cfg',
                            format='cfg')

        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        disp = 0.5 * np.sqrt(6.) / 3. * self._alat

        print(supercell_base[2, 2])

        #  crit1 = 20 * np.sqrt(3.) / 2. * 3.1648492
        crit1 = 40 * np.sqrt(3.) / 2. * 3.1648492

        count = 0
        for i in range(len(atom_position)):
            # and (atom_position[i, 2] < crit2):
            if (atom_position[i, 2] > crit1):
                atom_position[i, 0] += disp
                count += 1

        print(count)
        print(count / len(atom_position))

        atoms.set_positions(atom_position)
        self.write_lmp_config_data(atoms)

    def intro_kink_screw_dislocations(self, atoms, center1, center2, H=None, tilpnt=1. / 4., sign=1):
        # find how many layer of atoms #
        #  (layerList, distance) = \
        #  ase.utils.geometry.get_layers(atoms, [0, 0, 1], tolerance=0.01);
        boundary = 'ppp'

        lxid = 0  # glide plane
        lyid = 2  # y
        lzid = 1  # dislocation line

        cell = atoms.get_cell()
        nleft = tilpnt * cell[lzid, lzid]
        nrigh = (1 - tilpnt) * cell[lzid, lzid]

        if H is None:
            H = 0.5 * (center2[lxid] - center1[lxid])

        alpha = 0.1
        xc0, yc0 = center1[0], center1[1]

        # cut atoms (do not need for peierodic boundary condition ) ###
        if not boundary == 'ppp':
            # not using the ppp bounary condition ####
            # we have some other treaments ####
            lowx = 5 * (0.5 * (2.56656))
            higx = cell[0, 0] - lowx

            lowy = 5 * 2.22271
            higy = cell[1, 1] - lowy

            ##########################################################
            # the vol I cut here is
            # (lowx * y + (lowy * (x - 2 * lowx))) * 2 * z
            ##########################################################
            vol = np.linalg.det(cell)
            areaxy = (cell[0, 0] - 2 * lowx) * cell[1, 1]
            volcut = ((lowx * cell[1, 1] + lowy *
                       (cell[0, 0] - 2 * lowx)) * 2 * cell[2, 2])

            print(("vol[%f] - volcut[%f] = %f " % (vol, volcut, vol - volcut)))
            print(("area is [%f]" % (areaxy)))  # 29094.182232  A^2

            index_list = []
        cut = None  # only cut y direction ####
        #  shift   = True;
        shift = False

        # loop over atoms to add displacement #
        for atom in atoms:
            pos = atom.position

            if not boundary == 'ppp':
                # cut  10 A along x and y  #
                if cut == "xy":
                    if ((pos[lxid] < lowx) or (pos[lxid] > higx) or
                            (pos[lxid] > higy) or (pos[lxid] < lowy)):
                        index_list.append(atom.index)

                elif cut == "y":
                    if ((pos[lyid] > higy) or (pos[lyid] < lowy)):
                        index_list.append(atom.index)

            # calculate the h #####
            h = H * (1 + math.tanh(alpha * (pos[lzid] - nleft))) - \
                H * (1 + math.tanh(alpha * (pos[lzid] - nrigh)))

            xc = xc0 + h
            yc = yc0

            # calculate the displacement #####
            dx = pos[lxid] - xc
            dy = pos[lyid] - yc

            # add shift if peierodic along x   ######
            if shift is True:
                if pos[lxid] > xc:
                    # 0.3 or 0.2 for e1 = 1 -1 0
                    pos[lzid] += 0.25 * self.burger

            theta = np.arctan2(dy, dx)
            dz = sign * self.screw_coeff * theta
            pos[lzid] -= dz

            # assign the position ###
            atom.position = pos

        if not boundary == 'ppp':
            # delete boundary atoms ####
            del atoms[index_list]
        return atoms


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = bcc_kink()

    dispatcher = {'build': drv.construct_init_final,
                  'dipole': drv.cal_kink_pair_neb_pre,
                  'dipoleneb': cal_kink_pair_neb_dipole}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

    #  drv.cal_kink_pair_neb_dipole('initial')
    #  drv.cal_kink_pair_neb_dipole('final')
    # drv.cal_kink_pair_neb_dipole('each_interp')
