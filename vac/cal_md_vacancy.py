#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-11-23 09:46:18
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-02 21:00:24


from optparse import OptionParser
from glob import glob
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import shutil
import gn_config
import get_data
import gn_pbs
import atomman as am
import atomman.lammps as lmp
import md_pot_data


class cal_md_vacancy(gn_config.gnStructure, get_data.get_data, gn_pbs.gn_pbs):

    def __init__(self):
        self.pot = self.load_data("../BASICS/pot.dat")
        # self.pot = md_pot_data.va_pot.Nb_pbe
        # self.pot = md_pot_data.md_pot.Nb_eam 
        gn_config.gnStructure.__init__(self, self.pot)

    def bcc_vacancy_prep(self):
        atoms = self.set_bcc_convention().repeat((16, 16, 16))
        self.mymkdir("bulk")
        self.write_lmp_config_data(atoms, "bulk/lmp_init.txt")
        os.system("cp in.minimize bulk")
        atoms.pop(16 * 16 * 8)
        self.mymkdir("vacancy")
        self.write_lmp_config_data(atoms, "vacancy/lmp_init.txt")
        os.system("cp in.minimize vacancy")

    def bcc_vacancy_run(self):
        dirlist = ["bulk", "vacancy"]
        for locdir in dirlist:
            os.chdir(locdir)
            os.system("lmp_mpi -i in.minimize")
            os.chdir(os.pardir)

    def bcc_vacancy_cal(self):
        bulkeng = np.loadtxt("bulk/out")
        vacancyeng = np.loadtxt("vacancy/out")
        with open("bulk/lmp_init.txt", 'r') as fid:
            system = lmp.atom_data.load(data=fid.read(), units='metal')
        natoms = system.atoms.natoms
        print("vac formation energy",
              vacancyeng - (natoms - 1) * bulkeng / (natoms))

    def migration_prep(self):
        sz = 16
        # sz = 4
        perf_atoms = self.set_bcc_convention(size=[sz, sz, sz])
        vect = np.array([1, 1, 1]) * self.pot['lattice']

        pos1 = vect * sz * 0.5
        pos2 = vect * sz * 0.5 + 0.5 * vect
        pos3 = vect * 0.0

        dellist = []
        atoms = perf_atoms.copy()
        for atom in atoms:
            if np.isclose(pos1 - atom.position, np.array([0, 0, 0])).all() or \
                    np.isclose(pos2 - atom.position,
                               np.array([0, 0, 0])).all() or \
                    np.isclose(pos3 - atom.position,
                               np.array([0, 0, 0])).all():
                dellist.append(atom.index)
        print(dellist)
        del atoms[dellist]

        atoms1 = atoms.copy()
        atoms2 = atoms.copy()

        atoms1.append(ase.Atom('Nb', pos1))
        # atoms1.append(ase.Atom('Nb', pos3))
        # set this atom to be other type so we can fix it in lammps
        atoms1.append(ase.Atom('W', pos3))
        atoms2.append(ase.Atom('Nb', pos2))
        # atoms2.append(ase.Atom('Nb', pos3))
        atoms2.append(ase.Atom('W', pos3))

        # for vasp
        self.write_poscar_fix(atoms1, "posi")  # change it manually then
        self.write_poscar_fix(atoms2, "posf")

        # for md
        system, elements = am.convert.ase_Atoms.load(atoms1)
        lmp.atom_data.dump(system, "init.txt")
        shutil.copy("init.txt", "init0.txt")

        system, elements = am.convert.ase_Atoms.load(atoms2)
        lmp.atom_data.dump(system, "final.txt")
        shutil.copy("final.txt", "init1.txt")

    def migration_run(self):
        os.system("lmp_mpi -i in.init")
        os.system("lmp_mpi -i in.final")
        self.write_lmp_coords(ase.io.read(
            glob("bcc.final.*"), format='lammps-dump'), 'final.coord')
        os.system("rm bcc.final.*")

    def cal_vac_form(self):  # vasp PAW_PBE
        engy_vac = -2510.53396459
        num_atoms = 250
        engy_cohe = self._pot['ebcc']
        eform = engy_vac - engy_cohe * (num_atoms - 1)
        print(eform)

    def auto(self):
        drv.bcc_vacancy_prep()
        drv.bcc_vacancy_run()
        drv.bcc_vacancy_cal()


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_md_vacancy()
    dispatcher = {'fprep': drv.bcc_vacancy_prep,
                  'frun': drv.bcc_vacancy_run,
                  'fcal': drv.bcc_vacancy_cal,
                  'mprep': drv.migration_prep,
                  'auto': drv.auto}
    # 'cal': cal_vac_form
    # 'mrun': migration_run,
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
