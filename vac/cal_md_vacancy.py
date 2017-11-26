#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-11-23 09:46:18
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-23 09:46:46


from optparse import OptionParser
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import shutil
import glob
import md_pot_data
import gn_config
import get_data
import gn_lmp_infile
import gn_pbs

import atomman as am
import atomman.lammps as lmp


class cal_md_vacancy(gn_config.hcp,
                     gn_config.bcc,
                     gn_config.fcc,
                     get_data.get_data,
                     gn_pbs.gn_pbs,
                     gn_lmp_infile.gn_md_infile):

    def __init__(self):
        gn_lmp_infile.gn_md_infile.__init__(self)

        self._pot = md_pot_data.md_pot.Nb_eam
        self._element = self._pot['element']
        self._lat = self._pot['latbcc']
        self._pottype = self._pot['pair_style']
        self._potfile = self._pot['file']
        self.atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                                     [0, 1, 0],
                                                                     [0, 0, 1]],
                                                         latticeconstant=self._lat,
                                                         size=(1, 1, 1),
                                                         symbol=self._element,
                                                         pbc=(1, 1, 1))
        self.root = os.getcwd()
        return

    def bcc_vacancy_prep(self):
        atoms = self.atoms
        atoms = atoms.repeat((16, 16, 16))
        self.write_lmp_config_data(atoms, "bulk.txt")
        self.gn_md_input_vacancy("bulk.txt",
                                 self._pottype,
                                 self._potfile,
                                 self._element)

        self.mymkdir("bulk")
        os.system("mv bulk.txt bulk")
        os.system("cp dummy.lammps.ADP bulk")
        os.system("mv in.minimize bulk")

        atoms.pop(16 * 16 * 8)
        self.mymkdir("vacancy")
        self.write_lmp_config_data(atoms, "vacancy.txt")
        self.gn_md_input_vacancy("vacancy.txt",
                                 self._pottype,
                                 self._potfile,
                                 self._element)

        os.system("mv vacancy.txt vacancy")
        os.system("cp dummy.lammps.ADP vacancy")
        os.system("mv in.minimize vacancy")
        return

    def bcc_vacancy_run(self):
        dirlist = ["bulk", "vacancy"]
        for locdir in dirlist:
            os.chdir(locdir)
            os.system("lmp_mpi -i in.minimize")
            os.chdir(self.root)
        return

    def bcc_vacancy_cal(self):
        os.system("cp bulk/log.lammps .")
        bulkeng = self.md_get_final_energy()

        os.system("cp vacancy/log.lammps .")
        vacancyeng = self.md_get_final_energy()

        os.system("cp bulk/bulk.txt .")
        fid = open("bulk.txt", 'r')
        data = fid.read()

        fid.close()
        system = lmp.atom_data.load(data,
                                    units='metal',
                                    atom_style='atomic')
        natoms = system.atoms.natoms

        print "bulk", bulkeng
        print "vacancy", vacancyeng

        eng = vacancyeng - (natoms - 1) * bulkeng / (natoms)
        print eng
        return

    def migration_prep(self):
        perf_atoms = self.atoms
        sizen = 16
        perf_atoms = perf_atoms.repeat((sizen, sizen, sizen))
        alat = self._lat
        vect = np.array([1, 1, 1]) * alat

        pos1 = vect * sizen * 0.5
        pos2 = vect * sizen * 0.5 + 0.5 * vect
        pos3 = vect * 0.0
        dellist = []
        atoms = perf_atoms.copy()
        for atom in atoms:
            if ((pos1 - atom.position) == ([0, 0, 0])).all() or \
                    ((pos2 - atom.position) == ([0, 0, 0])).all() or \
                    ((pos3 - atom.position) == ([0, 0, 0])).all():
                dellist.append(atom.index)
        print dellist
        del atoms[dellist]

        atoms1 = atoms.copy()
        atoms2 = atoms.copy()

        atoms1.append(ase.Atom('Nb', pos1))
        atoms1.append(ase.Atom('W', pos3))
        atoms2.append(ase.Atom('Nb', pos2))
        atoms2.append(ase.Atom('W', pos3))

        system, elements = am.convert.ase_Atoms.load(atoms1)
        lmp.atom_data.dump(system, "init.txt")
        shutil.copy("init.txt", "init0.txt")

        system, elements = am.convert.ase_Atoms.load(atoms2)
        lmp.atom_data.dump(system, "final.txt")
        shutil.copy("final.txt", "init1.txt")
        return

    def migration_run(self):
        os.system("lmp_mpi -i in.init")
        os.system("lmp_mpi -i in.final")
        filein = glob.glob("bcc.final.*")
        ase_atoms = ase.io.read(filein[-1],
                                format='lammps-dump')
        self.write_lmp_coords(ase_atoms,
                              filename='final.coord')
        os.system("rm bcc.final.*")
        return

    def cal_vac_form(self):
        # vasp PAW_PBE
        engy_vac = -2510.53396459
        num_atoms = 250
        engy_cohe = self._pot['ebcc']
        eform = engy_vac - engy_cohe * (num_atoms - 1)
        print eform
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_md_vacancy()

    if options.mtype.lower() == 'fprep':
        drv.bcc_vacancy_prep()

    if options.mtype.lower() == 'frun':
        drv.bcc_vacancy_run()

    if options.mtype.lower() == 'fcal':
        drv.bcc_vacancy_cal()

    if options.mtype.lower() == 'fauto':
        drv.bcc_vacancy_prep()
        drv.bcc_vacancy_run()
        drv.bcc_vacancy_cal()

    if options.mtype.lower() == 'mprep':
        drv.migration_prep()

    if options.mtype.lower() == 'mrun':
        drv.migration_run()

    if options.mtype.lower() == 'cal':
        drv.cal_vac_form()
