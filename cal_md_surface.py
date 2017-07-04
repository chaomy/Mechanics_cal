#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_md_surface.py
#
###################################################################
#
# Purpose :  multithreads to calculate md surface
#
# Creation Date :
# Last Modified : Thu Mar 30 23:56:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
from multiprocessing import Pool

try:
    import gn_config
    import get_data
    import gn_lmp_infile
    import output_data
    import md_pot_data

except ImportError:
    print "error during import"


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_surface.run_lmp_minimize(*arg, **kwarg)


class cal_md_surface(gn_config.bcc,
                     gn_config.fcc,
                     gn_config.hcp,
                     get_data.get_data,
                     gn_lmp_infile.gn_md_infile,
                     output_data.output_data):

    def __init__(self):
        self.pot = md_pot_data.md_pot.Nb_adp
        get_data.get_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        output_data.output_data.__init__(self)

        self._surface_type = '100'
        self._surface_element = self.pot['element']
        self._surface_lattice_constant = self.pot['lattice']
        self._structure = self.pot['structure']
        self._surface_potential = self.pot['file']

        if self._structure == 'bcc':
            gn_config.bcc.__init__(self, self.pot)
        elif self._structure == 'fcc':
            gn_config.fcc.__init__(self, self.pot)
        elif self._structure == 'hcp':
            gn_config.hcp.__init__(self, self.pot)

        self.set_lattce_constant(self._surface_lattice_constant)
        self.set_element(self._surface_element)
        self.set_config_file_format("lmp")
        self.config_file = "lmp_init.txt"

        self.root_dir = os.getcwd()
        return

    def set_md__surface_potential(self, in_potential):
        self._surface_potential = in_potential
        return

    def set_surface_type(self, surface_type):
        self._surface_type = surface_type
        return

    def gn_surface_atoms(self):
        if self._surface_type == '100':
            self._surface_direction = [[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 20))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._surface_type == '110':
            self._surface_direction = [[-1, 1,  0],
                                       [0,  0,  1],
                                       [1,  1,  0]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._surface_type == '111':
            self._surface_direction = [[1,  1, -2],
                                       [-1, 1,  0],
                                       [1,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._surface_type == '112':
            self._surface_direction = [[1,  1, -2],
                                       [-1, 1,  0],
                                       [1,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

    def gn_bulk_atoms(self):
        if self._surface_type == '100':
            self._surface_direction = [[1,  0,  0],
                                       [0,  1,  0],
                                       [0,  0,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            return atoms

        elif self._surface_type == '110':
            self._surface_direction = [[-1, 1,  0],
                                       [0,  0,  1],
                                       [1,  1,  0]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 11))
            return atoms

        elif self._surface_type == '111':
            self._surface_direction = [[1, 1, -2],
                                       [-1, 1, 0],
                                       [1, 1, 1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 10))
            return atoms

        elif self._surface_type == '112':
            self._surface_direction = [[1, 1, -2],
                                       [-1, 1, 0],
                                       [1, 1, 1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 10))
            return atoms

    def run_lmp_minimize(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.minimize")
        os.chdir(self.root_dir)
        return

    def multi_thread_surface(self):
        dir_list = []
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            dir_list.append('dir-surf-%s' % (loop_list[i]))
            dir_list.append('dir-bulk-%s' % (loop_list[i]))
        num_threads = len(dir_list)
        pool = Pool(processes=num_threads)
        pool.map(unwrap_self_run_lammps,
                 zip([self] * num_threads,
                     dir_list))
        return

    def loop_cal_surface_data(self):
        if os.path.isfile("surface_energy.dat"):
            os.system("mv surface_energy.dat surface_backup.dat")
        loop_list = ['100', '110', '111']
        surface_energy_list = []
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            surface_energy_list.append(self.cal_surface_energy())
        return surface_energy_list

    def loop_prepare_suface(self):
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            self.prepare_md_surface()
            self.prepare_md_bulk()
        return

    def prepare_md_surface(self):
        dir_surf = 'dir-surf-%s' % (self._surface_type)

        if not os.path.isdir(dir_surf):
            os.mkdir(dir_surf)
        os.chdir(dir_surf)

        atoms = self.gn_surface_atoms()

        os.system("cp  ../%s  ." % (self._surface_potential))
        self.write_lmp_config_data(atoms)
        self.gn_md_minimize("lmp_init.txt",
                            self._surface_potential,
                            self._surface_element)

        os.chdir(self.root_dir)
        return

    def prepare_md_bulk(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)

        if not os.path.isdir(dir_bulk):
            os.mkdir(dir_bulk)
        os.chdir(dir_bulk)

        atoms = self.gn_bulk_atoms()

        os.system("cp  ../%s  ." % (self._surface_potential))
        self.write_lmp_config_data(atoms)
        self.gn_md_minimize("lmp_init.txt",
                            self._surface_potential,
                            self._surface_element)

        os.chdir(self.root_dir)
        return

    def cal_surface_energy(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        dir_surf = 'dir-surf-%s' % (self._surface_type)

        os.chdir(dir_bulk)
        energy_b = self.md_get_final_energy()
        super_cell = self.md_get_cell()
        xy_area = self.cal_poscar_xy_area(super_cell)
        os.chdir(self.root_dir)

        os.chdir(dir_surf)
        energy_s = self.md_get_final_energy()
        os.chdir(self.root_dir)

        surface_e = 0.5 * (energy_s - energy_b) / xy_area

        with open("surface_energy.dat", 'a') as fid:
            fid.write("surface energy %s %6.4f ev/A \n"
                      % (self._surface_type, surface_e))
            fid.write("surface energy %s %6.4f J/m \n"
                      % (self._surface_type, surface_e * self.ev_angstrom2_to_j_m2))
        return (surface_e * self.ev_angstrom2_to_j_m2)


if __name__ == '__main__':
    Job = cal_md_surface()

    Job.loop_prepare_suface()
    Job.multi_thread_surface()
    print Job.loop_cal_surface_data()