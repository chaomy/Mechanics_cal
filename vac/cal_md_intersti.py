#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-28 16:31:13


from optparse import OptionParser
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import md_pot_data
import atomman as am
import atomman.lammps as lmp
import atomman.unitconvert as uc
import gn_config
import get_data
import gn_lmp_infile
import gn_pbs


class cal_md_intersti(gn_config.gnStructure,
                      get_data.get_data,
                      gn_pbs.gn_pbs,
                      gn_lmp_infile.gn_md_infile):

    def __init__(self):
        gn_lmp_infile.gn_md_infile.__init__(self)
        self.pot = self.load_data("../BASICS/pot.dat")
        gn_config.gnStructure.__init__(self, self.pot)
        self._unitatoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                                          [0, 1, 0],
                                                                          [0, 0, 1]],
                                                              latticeconstant=self.pot[
                                                                  'lattice'],
                                                              size=(1, 1, 1),
                                                              symbol=self.pot[
                                                                  "element"],
                                                              pbc=(1, 1, 1))

    def cal_dumbbell_100(self, atoms, dirname="dumb100"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        #  b = 0.25
        x = (np.sqrt(10) - 1) / 3.
        b = 0.5 * (1 - x)

        # add atoms #
        pos1 = np.array([0.5, 0.5, b]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))

        pos2 = np.array([0.5, 0.5, 1 - b]) * alat
        atoms.append(ase.atom.Atom('W', position=pos2))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("energy increment :", engyincre)
        # return engyincre

    def cal_dumbbell_110(self, atoms, dirname="dumb110"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        b1 = (4 * np.sqrt(2) - np.sqrt(11)) / (6 * np.sqrt(2))
        b2 = 1 - b1

        # add atoms #
        pos1 = np.array([b1, b1, 0.5]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))
        pos2 = np.array([b2, b2, 0.5]) * alat
        atoms.append(ase.atom.Atom('W', position=pos2))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("110 dumbel interstitial: ", engyincre)
        # return engyincre

    def cal_dumbbell_111(self, atoms, dirname="dumb111"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        b1 = 1. / 3.
        b2 = 2. / 3.

        pos1 = np.array([b1, b1, b1]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))

        pos2 = np.array([b2, b2, b2]) * alat
        atoms.append(ase.atom.Atom('W', position=pos2))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("111 dumbel interstitial: ", engyincre)
        # return engyincre

    def cal_crowdion(self, atoms, dirname="crowdion"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        b1 = 0.25

        pos1 = np.array([b1, b1, b1]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("energy crowdion: ", engyincre)
        # return engyincre

    def cal_octahedral(self, atoms, dirname="octahedral"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        b = 0.5
        pos1 = np.array([b, b, 0]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("energy octahedral: ", engyincre)
        # return engyincre

    def cal_tetrahedral(self, atoms, dirname="eoctahedral"):
        alat = self.pot['lattice']
        eperatom = self.pot['ebcc']

        pos1 = np.array([0.25, 0.5, 0]) * alat
        atoms.append(ase.atom.Atom('W', position=pos1))

        self.mymkdir(dirname)
        self.write_lmp_config_data(atoms, "{}/init.txt".format(dirname))

        # engyinter = self.write_and_run(dirname, atoms)
        # engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        # print("energy octahedral: ", engyincre)
        # return engyincre

    def prep_interstitials(self):  # reported size 25 x 25 x 25
        alat = self.pot['lattice']
        atoms = self._unitatoms.copy()
        sizen = 25
        atoms = atoms.repeat((sizen, sizen, sizen))
        pos = np.array([alat,  alat,  alat]) * 0.5

        for atom in atoms:
            if ((pos - atom.position) == ([0, 0, 0])).all():
                index = atom.index
                print("del {}  atom".format(index))
        del atoms[index]

        atomsper = self._unitatoms.copy().repeat((sizen, sizen, sizen))

        self.cal_dumbbell_100(atoms.copy())
        self.cal_dumbbell_110(atoms.copy())
        self.cal_dumbbell_111(atoms.copy())

        self.cal_crowdion(atomsper.copy())
        self.cal_octahedral(atomsper.copy())
        self.cal_tetrahedral(atomsper.copy())

    def run(self):
        dls = ["dumb100", "dumb110", "dumb111",
               "crowdion", "eoctahedral", "octahedral"]
        for dd in dls:
            os.chdir(dd)
            os.system("cp ../in.init .")
            os.system("mpirun lmp_mpi -i in.init")
            os.chdir(os.pardir)

    def load(self):
        dls = ["dumb100", "dumb110", "dumb111",
               "crowdion", "eoctahedral", "octahedral"]
        for dd in dls:
            data = np.loadtxt("{}/out".format(dd))
            print(dd, data[0] - data[1] * self.pot["ebcc"])

        va_pbe = {"100": -2521.89426063,
                  "111": -2522.80887145}


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()

    drv = cal_md_intersti()

    dispatcher = {'prep': drv.prep_interstitials,
                  'run': drv.run,
                  'load': drv.load}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
