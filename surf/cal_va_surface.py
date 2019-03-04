#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2019-03-03 21:57:02


from optparse import OptionParser
import os
import ase.io


class cal_va_surface(object):

    def prep_va_inputs(self, mdir):
        # INCAR #
        self.set_incar_type('ab_md')
        self.set_accuracy(1e-5)
        self.write_incar()
        # KPOINTS #
        self.set_intype('gamma')
        self.write_kpoints()
        # va.pbs #
        self.set_pbs(mdir, 'va')
        return

    def prep_va_surface(self, dirtag='dir', mtype='relax'):
        configs = ['surf', 'bulk']
        atomsl = self.gn_surface_atoms()
        for config, atoms in zip(configs, atomsl):
            if self.mgsf in ['x100z100']:
                atoms = atoms.repeat((2, 2, 1))
            if self.mgsf in ['x110z110']:
                atoms = atoms.repeat((2, 2, 1))
            if self.mgsf in ['x112z111']:
                atoms = atoms.repeat((1, 2, 1))
            mdir = '{}-{}-{}'.format(dirtag, self.mgsf, config)
            self.mymkdir(mdir)
            os.chdir(mdir)
            ase.io.write('POSCAR', images=atoms, format='vasp')
            self.prep_va_inputs(mdir)
            os.system("cp ../POTCAR  .")
            os.chdir(os.pardir)
        return

    def prepare_vasp_bulk(self):
        dir_bulk = 'dir-bulk-%s' % (self._surface_type)
        self.set_diff_kpoints(self._surface_kpoints)
        os.chdir(os.pardir)
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
        mdir = 'dir-surf-%s' % (self._surface_type)

        os.chdir(dir_bulk)
        (energy_b, stress, volume) = self.vasp_energy_stress_vol()
        (atom_number, supercell_base, comment, atom_position) = \
            self.read_vasp_poscar()
        os.chdir(os.pardir)

        os.chdir(mdir)
        (energy_s, stress, volume) = self.vasp_energy_stress_vol()
        os.chdir(os.pardir)

        xy_area = self.cal_poscar_xy_area(supercell_base)
        surface_e = 0.5 * (energy_s - energy_b) / xy_area

        with open("surface_energy.dat", 'a') as fid:
            fid.write("surface energy %s %6.4f ev/A \n"
                      % (self._surface_type, surface_e))
            fid.write("surface energy %s %6.4f J/m \n"
                      % (self._surface_type, surface_e * self.ev_angstrom_to_j_m))
        return

    def loop_sub_jobs(self):
        dir_list = glob.glob("dir[_-]*")
        for i in range(len(dir_list)):
            os.chdir(dir_list[i])
            os.system("qsub va.pbs")
            os.chdir(os.pardir)


if __name__ == '__main__':
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
    Job = cal_surface(in_structure='bcc',
                      surface_type='100',
                      lattice_constant=3.322404,
                      in_element='Nb',
                      in_kpoints=[18, 18, 1])

    if options.mtype == 'prep':
        Job.loop_prepare_suface_run()

    if options.mtype == 'clc':
        Job.loop_cal_surface_data()
