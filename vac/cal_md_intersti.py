#!/usr/bin/env python
# encoding: utf-8


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


class cal_md_intersti(gn_config.hcp,
                      gn_config.bcc,
                      gn_config.fcc,
                      get_data.get_data,
                      gn_pbs.gn_pbs,
                      gn_lmp_infile.gn_md_infile):

    def __init__(self):
        gn_lmp_infile.gn_md_infile.__init__(self)
        self.pot = md_pot_data.md_pot.Nb_adp

        self._lat = 3.30789893314849
        self._pottype = 'adp'
        self._coherent_eam = -7.09123298005381  # ev / atom

        self._element = self.pot['element']
        self._lat = self.pot['latbcc']
        self._pottype = self.pot['pair_style']
        self._potfile = self.pot['file']

        self._unitatoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                                          [0, 1, 0],
                                                                          [0, 0, 1]],
                                                              latticeconstant=self._lat,
                                                              size=(1, 1, 1),
                                                              symbol=self._element,
                                                              pbc=(1, 1, 1))
        self._root = os.getcwd()
        self._run_lmp = "lmp_mpi -i in.init"

    def write_and_run(self, dirname, atoms):
        system, elements = am.convert.ase_Atoms.load(atoms)
        lmp.atom_data.dump(system, "init.txt")

        self.mymkdir(dirname)
        os.system("mv init.txt  {}".format(dirname))
        os.system("cp in.init {}".format(dirname))

        os.chdir(dirname)
        os.system(self._run_lmp)
        engyinter = self.md_get_final_energy()
        os.chdir(self._root)
        return engyinter

    def cal_dumbbell_100(self, atoms):
        dirname = "dumb100"
        alat = self._lat
        eperatom = self._coherent_eam

        #  b = 0.25
        x = (np.sqrt(10) - 1) / 3.
        b = 0.5 * (1 - x)

        # add atoms #
        pos1 = np.array([0.5, 0.5, b]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        pos2 = np.array([0.5, 0.5, 1 - b]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos2))

        #  write and run #
        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("energy increment :", engyincre)
        return engyincre

    def cal_dumbbell_110(self, atoms):
        dirname = "dumb110"
        alat = self._lat
        eperatom = self._coherent_eam

        b1 = (4 * np.sqrt(2) - np.sqrt(11)) / (6 * np.sqrt(2))
        b2 = 1 - b1

        # add atoms #
        pos1 = np.array([b1, b1, 0.5]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        pos2 = np.array([b2, b2, 0.5]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos2))

        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("110 dumbel interstitial: ", engyincre)
        return engyincre

    def cal_dumbbell_111(self, atoms):
        dirname = "dumb111"
        alat = self._lat
        eperatom = self._coherent_eam

        b1 = 1. / 3.
        b2 = 2. / 3.

        # add atoms #
        pos1 = np.array([b1, b1, b1]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        pos2 = np.array([b2, b2, b2]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos2))

        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("111 dumbel interstitial: ", engyincre)
        return engyincre

    def cal_crowdion(self, atoms):
        dirname = "crowdion"
        alat = self._lat
        eperatom = self._coherent_eam

        b1 = 0.25
        # add atoms #
        pos1 = np.array([b1, b1, b1]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("energy crowdion: ", engyincre)
        return engyincre

    def cal_octahedral(self, atoms):
        dirname = "octahedral"
        alat = self._lat
        eperatom = self._coherent_eam

        b = 0.5
        # add atoms #
        pos1 = np.array([b, b, 0]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("energy octahedral: ", engyincre)
        return engyincre

    def cal_tetrahedral(self, atoms):
        dirname = "eoctahedral"
        alat = self._lat
        eperatom = self._coherent_eam

        pos1 = np.array([0.25, 0.5, 0]) * alat
        atoms.append(ase.atom.Atom('W',
                                   position=pos1))

        engyinter = self.write_and_run(dirname, atoms)
        engyincre = engyinter - eperatom * (atoms.get_number_of_atoms())
        print("energy octahedral: ", engyincre)
        return engyincre

    def prep_interstitials(self):
        # reported size 25 x 25 x 25
        alat = self._lat
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

        edumbb100 = self.cal_dumbbell_100(atoms.copy())
        edumbb110 = self.cal_dumbbell_110(atoms.copy())
        edumbb111 = self.cal_dumbbell_111(atoms.copy())

        ecrowdion = self.cal_crowdion(atomsper.copy())
        eoctahedral = self.cal_octahedral(atomsper.copy())
        etetrahedral = self.cal_tetrahedral(atomsper.copy())
        print("100 110 111", edumbb100, edumbb110, edumbb111)
        print("ecrowdion octahedral etetrahedral",
              ecrowdion, eoctahedral, etetrahedral)
        return

    def loop_collect_r(self):  # total is 14 ####
        tag = "glob"
        if tag == "given":
            for i in range(6):
                r_min = 2.0 + 0.2 * i
                r_max = 2.2 + 0.2 * i
                #  self.generate_rand_given_distance(r_min, r_max);
                dirname = "r_%.2f_%.2f" % (r_min, r_max)
                os.chdir(dirname)
                os.system("cp OUTCAR OUTCAR_ref")
                os.system("vasp2force -f OUTCAR_ref >> ../dummy.config_r")
                os.chdir(self.root)

        elif tag == "glob":
            file_list = glob.glob("r_*")
            for i in range(len(file_list)):
                dirname = file_list[i]
                os.chdir(dirname)
                os.system("cp OUTCAR OUTCAR_ref")
                os.system("vasp2force -f OUTCAR_ref >> ../dummy.config_r")
                os.chdir(self.root)
        return

    def loop_prep_r(self):
        #### total is 14 ####
        for i in range(0, 2):
            r_min = 2.0 + 0.2 * i
            r_max = r_min + 0.2

            dirname = "r_%.2f_%.2f" % (r_min, r_max)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)

            os.chdir(dirname)
            self.generate_rand_given_distance(r_min, r_max, dirname)

            #  os.system("qsub va.pbs");
            os.chdir(self.root)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_md_intersti()

    if options.mtype.lower() == 'prep':
        drv.prep_interstitials()
