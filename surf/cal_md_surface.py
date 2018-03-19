#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-18 09:45:55


from multiprocessing import Pool
from optparse import OptionParser
import os
import gn_config
import get_data
import gsf_data
import gn_lmp_infile
import md_pot_data
import output_data


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_surface.run_lmp_minimize(*arg, **kwarg)


class cal_md_surface(gn_config.bcc,
                     get_data.get_data,
                     gn_lmp_infile.gn_md_infile,
                     output_data.output_data):

    def __init__(self, pot=None):
        if pot is None:
            self.pot = md_pot_data.md_pot.Nb_meamc
        else:
            self.pot = pot
        get_data.get_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        output_data.output_data.__init__(self)
        self._surface_type = '100'
        gn_config.bcc.__init__(self, self.pot)
        self.set_config_file_format("lmp")
        self.config_file = "lmp_init.txt"

    def set_surface_type(self, surface_type):
        self._surface_type = surface_type

    def gn_surface_atoms(self):
        if self._surface_type in ['100']:
            tag = 'x100z100'
        elif self._surface_type in ['110']:
            tag = 'x110z110'
        elif self._surface_type in ['111']:
            tag = 'x112z111'
        atoms = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[tag],
            in_size=gsf_data.gsfsize[tag])
        for i in range(gsf_data.gsfpopn[tag]):
            atoms.pop()
        return atoms

    def gn_bulk_atoms(self):
        if self._surface_type in ['100']:
            tag = 'x100z100'
        elif self._surface_type in ['110']:
            tag = 'x110z110'
        elif self._surface_type in ['111']:
            tag = 'x112z111'
        atoms = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[tag],
            in_size=gsf_data.bulksize[tag])
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

    def prepare_md_surface(self):
        dir_surf = 'dir-surf-%s' % (self._surface_type)
        self.mymkdir(dir_surf)
        self.write_lmp_config_data(self.gn_surface_atoms())
        os.system("cp lmp_init.txt {}".format(dir_surf))
        os.system('cp in.minimize  {}'.format(dir_surf))

    def prepare_md_bulk(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        self.mymkdir(dir_bulk)
        self.write_lmp_config_data(self.gn_bulk_atoms())
        os.system("cp lmp_init.txt {}".format(dir_bulk))
        os.system('cp in.minimize  {}'.format(dir_bulk))

    def cal_surface_energy(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        dir_surf = 'dir-surf-%s' % (self._surface_type)

        os.chdir(dir_bulk)
        energy_b = self.md_get_final_energy()
        super_cell = self.md_get_cell()
        xy_area = self.cal_poscar_xy_area(super_cell)
        os.chdir(os.pardir)

        os.chdir(dir_surf)
        energy_s = self.md_get_final_energy()
        os.chdir(os.pardir)

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
