#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-06 14:07:24

from multiprocessing import Pool
from optparse import OptionParser
import os
import gn_config
import get_data
import gn_lmp_infile
import md_pot_data
import output_data


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_surface.run_lmp_minimize(*arg, **kwarg)


class cal_md_surface(gn_config.bcc,
                     get_data.get_data,
                     gn_lmp_infile.gn_md_infile,
                     output_data.output_data):

    def __init__(self):
        # self.pot = md_pot_data.md_pot.Nb_eam
        self.pot = self.load_data('pot.dat')
        get_data.get_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        output_data.output_data.__init__(self)

        self._surface_type = '100'
        gn_config.bcc.__init__(self, self.pot)
        self.set_lattce_constant(self.pot['latbcc'])
        self.set_config_file_format("lmp")
        self.config_file = "lmp_init.txt"
        self.root_dir = os.getcwd()
        return

    def set_surface_type(self, surface_type):
        self._surface_type = surface_type
        return

    def gn_surface_atoms(self):
        if self._surface_type == '100':
            self._surdir = [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 20))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._surface_type == '110':
            self._surdir = [[-1, 1, 0],
                            [0, 0, 1],
                            [1, 1, 0]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 14))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._surface_type == '111':
            self._surdir = [[1, 1, -2],
                            [-1, 1, 0],
                            [1, 1, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._surface_type == '112':
            self._surdir = [[1, 1, -2],
                            [-1, 1, 0],
                            [1, 1, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

    def gn_bulk_atoms(self):
        if self._surface_type == '100':
            self._surdir = [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 14))
            return atoms

        elif self._surface_type == '110':
            self._surdir = [[-1, 1, 0],
                            [0, 0, 1],
                            [1, 1, 0]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 11))
            return atoms

        elif self._surface_type == '111':
            self._surdir = [[1, 1, -2],
                            [-1, 1, 0],
                            [1, 1, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 10))
            return atoms

        elif self._surface_type == '112':
            self._surdir = [[1, 1, -2],
                            [-1, 1, 0],
                            [1, 1, 1]]
            atoms = self.set_bcc_convention(
                in_direction=self._surdir,
                in_size=(1, 1, 10))
            return atoms

    def run_lmp_minimize(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.minimize")
        os.chdir(os.pardir)
        return

    def loop_run_surface(self):
        dir_list = []
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            dir_list.append('dir-surf-%s' % (loop_list[i]))
            dir_list.append('dir-bulk-%s' % (loop_list[i]))
        for mdir in dir_list:
            self.run_lmp_minimize(mdir)
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
                 zip([self] * num_threads, dir_list))
        return

    def loop_cal_surface_data(self):
        loop_list = ['100', '110', '111']
        surface_energy_list = []
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            surface_energy_list.append(self.cal_surface_energy())
        return surface_energy_list

    def loop_prepare_surface(self):
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            self.prepare_md_surface()
            self.prepare_md_bulk()
        return

    def prepare_md_surface(self):
        dir_surf = 'dir-surf-%s' % (self._surface_type)
        self.mymkdir(dir_surf)
        os.chdir(dir_surf)
        atoms = self.gn_surface_atoms()
        os.system("cp  ../%s  ." % (self.pot['file']))
        self.write_lmp_config_data(atoms)
        os.system('cp  ../in.minimize  .')
        os.chdir(self.root_dir)
        return

    def prepare_md_bulk(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        self.mymkdir(dir_bulk)
        os.chdir(dir_bulk)
        atoms = self.gn_bulk_atoms()
        os.system("cp  ../%s  ." % (self.pot['file']))
        os.system('cp  ../in.minimize  .')
        self.write_lmp_config_data(atoms)
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
                      % (self._surface_type,
                         surface_e * self.ev_angstrom2_to_j_m2))
        return (surface_e * self.ev_angstrom2_to_j_m2)

    def wrap_run(self):
        self.loop_prepare_surface()
        self.multi_thread_surface()
        print self.loop_cal_surface_data()
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_md_surface()
    dispatcher = {'run': drv.wrap_run}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
