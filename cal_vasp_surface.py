#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_vasp_surface.py
#
###################################################################
#
# Purpose : calculate surface energy by vasp
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
import glob
from optparse import OptionParser

try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs

except ImportError:
    print "error during import"


class cal_surface(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  get_data.get_data,
                  gn_kpoints.gn_kpoints,
                  gn_incar.gn_incar,
                  gn_pbs.gn_pbs):

    def __init__(self,
                 in_structure='bcc',
                 surface_type='100',
                 lattice_constant=3.17,
                 in_element='W',
                 in_kpoints=[18, 18, 1]):

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)

        self._surface_type = surface_type
        self.surface_type = '100'  # default
        self._surface_element = in_element
        self._surface_lattice_constant = lattice_constant
        self._surface_kpoints = in_kpoints
        self._structure = in_structure

        if self._structure == 'bcc':
            gn_config.bcc.__init__(self,
                                   self._surface_element,
                                   self._surface_lattice_constant)

        elif self._structure == 'fcc':
            gn_config.fcc.__init__(self,
                                   self._surface_element,
                                   self._surface_lattice_constant)

        elif self._structure == 'hcp':
            gn_config.hcp.__init__(self,
                                   self._surface_element,
                                   self._surface_lattice_constant)

        self.exe = "mpirun vasp"

        self.root_dir = os.getcwd()
        self.set_accuracy(1e-5)
        return

    def set_surface_type(self, surface_type):
        self._surface_type = surface_type
        return

    def set_surface_norm_direction(self):
        if self._surface_type == '100':
            self._surface_direction = [[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]

        elif self._surface_type == '110':
            self._surface_direction = [[-1, 1, 0],
                                       [0,  0, 1],
                                       [1,  1, 0]]

        elif self._surface_type == '111':
            self._surface_direction = [[1,  1, -2],
                                       [-1, 1, 0],
                                       [1,  1, 1]]
        return

    def gn_surface_100(self):
        self.set_surface_type('100')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))  # default 18
        for i in range(8):
            atoms.pop()
        return atoms

    def gn_surface_110(self):
        self.set_surface_type('110')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 12))
        for i in range(12):
            atoms.pop()
        return atoms

    def gn_surface_111(self):
        self.set_surface_type('111')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 16))
        for i in range(24):
            atoms.pop()
        return atoms

    def gn_bulk_100(self):
        self.set_surface_type('100')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))
        return atoms

    def gn_bulk_110(self):
        self.set_surface_type('110')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 11))
        return atoms

    def gn_bulk_111(self):
        self.set_surface_type('111')
        self.set_surface_norm_direction()
        atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 10))
        return atoms

    def get_potcar(self):
        os.system("cp ../../POTCAR .")
        return

    def prepare_vasp_inputs(self, dir_name):
        #### INCAR ####
        self.set_incar_type('dft')
        self.set_accuracy(1e-5)
        self.write_incar()

        #### KPOINTS ####
        self.set_intype('gamma')
        self.write_kpoints()

        #### va.pbs ####
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_pbs_type('va')
        self.set_wall_time(120)
        self.set_job_title(dir_name)
        self.set_main_job("mpirun vasp")
        self.write_pbs()

        #  self.get_potcar()
        return

    def prepare_vasp_surface(self):
        dir_surf = 'dir-surf-%s' % (self._surface_type)
        self.set_diff_kpoints(self._surface_kpoints)

        if not os.path.isdir(dir_surf):
            os.mkdir(dir_surf)
        os.chdir(dir_surf)

        if self._surface_type == '100':
            atoms = self.gn_surface_100()
        elif self._surface_type == '110':
            atoms = self.gn_surface_110()
        elif self._surface_type == '111':
            atoms = self.gn_surface_111()
        self.write_poscar(atoms)
        self.prepare_vasp_inputs(dir_surf)

        os.chdir(self.root_dir)
        return

    def prepare_vasp_bulk(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        self.set_diff_kpoints(self._surface_kpoints)

        if not os.path.isdir(dir_bulk):
            os.mkdir(dir_bulk)
        os.chdir(dir_bulk)

        if self._surface_type == '100':
            atoms = self.gn_bulk_100()
        elif self._surface_type == '110':
            atoms = self.gn_bulk_110()
        elif self._surface_type == '111':
            atoms = self.gn_bulk_111()

        self.set_incar_type('dft')
        self.write_incar()
        self.write_kpoints()
        self.write_poscar(atoms)

        self.set_pbs_type('va')
        self.set_wall_time(120)
        self.set_job_title(dir_bulk)
        self.set_main_job("mpirun vasp")
        self.write_pbs()
        #  self.get_potcar()

        os.chdir(self.root_dir)
        return

    def loop_prepare_suface_run(self):
        #  loop_list = ['100', '110', '111']
        loop_list = ['110']
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            self.prepare_vasp_surface()
            self.prepare_vasp_bulk()
        return

    def loop_cal_surface_data(self):
        if os.path.isfile("surface_energy.dat"):
            os.system("mv surface_energy.dat surface_backup.dat")
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            loc_type = loop_list[i]
            self.set_surface_type(loc_type)
            self.cal_surface_energy()
        return

    def cal_surface_energy(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        dir_surf = 'dir-surf-%s' % (self._surface_type)

        os.chdir(dir_bulk)
        (energy_b, stress, volume) = self.vasp_energy_stress_vol()
        (atom_number, supercell_base, comment, atom_position) = \
            self.read_vasp_poscar()
        os.chdir(self.root_dir)

        os.chdir(dir_surf)
        (energy_s, stress, volume) = self.vasp_energy_stress_vol()
        os.chdir(self.root_dir)

        xy_area = self.cal_poscar_xy_area(supercell_base)
        surface_e = 0.5 * (energy_s - energy_b) / xy_area

        with open("surface_energy.dat", 'a') as fid:
            fid.write("surface energy %s %6.4f ev/A \n"
                      % (self._surface_type, surface_e))
            fid.write("surface energy %s %6.4f J/m \n"
                      % (self._surface_type, surface_e * self.ev_angstrom_to_j_m))
        return

    def loop_sub_jobs(self):
        dir_list = glob.glob("dir-*")
        for i in range(len(dir_list)):
            os.chdir(dir_list[i])
            os.system("qsub va.pbs")
            os.chdir(self.root_dir)
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
    ##### W  3.1711    #####
    ##### Nb 3.322404  #####
    Job = cal_surface(in_structure='bcc',
                      surface_type='100',
                      lattice_constant=3.322404,
                      in_element='Nb',
                      in_kpoints=[18, 18, 1])

    if options.mtype == 'prep':
        Job.loop_prepare_suface_run()

    if options.mtype == 'sub':
        Job.loop_sub_jobs()

    if options.mtype == 'clc':
        Job.loop_cal_surface_data()
