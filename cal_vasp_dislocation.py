#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_vasp_dislocation.py
#
###################################################################
#
# Purpose :
#
# Creation Date : Tue Apr 11 15:36:33 2017
# Last Modified : Sat Apr  1 23:15:41 2017
# Created By    : Chaoming Yang
#
###################################################################


import glob
import os
import numpy as np
import shutil
from optparse import OptionParser

try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import Intro_vasp
    import gn_pbs
    import dd_map

except ImportError:
    print "error during import"


class vasp_dislocation(gn_config.bcc,
                       gn_config.fcc,
                       gn_config.hcp,
                       get_data.get_data,
                       gn_kpoints.gn_kpoints,
                       gn_incar.gn_incar,
                       gn_pbs.gn_pbs,
                       Intro_vasp.vasp_change_box):

    def __init__(self, structure='bcc'):
        gn_incar.gn_incar.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        gn_pbs.gn_pbs.__init__(self)

        self._lattice_constant = 3.1711  # W dft
        self.lattice_constant = 3.322404  # Nb dft
        #  self._lattice_constant = 3.143390   # w_eam4

        Intro_vasp.vasp_change_box.__init__(self, self._lattice_constant)
        self.structure = structure

        if self.structure == 'bcc':
            gn_config.bcc.__init__(self)
        elif self.structure == 'fcc':
            gn_config.fcc.__init__(self)
        elif self.structure == 'hcp':
            gn_config.hcp.__init__(self)

        self.root_dir = os.getcwd()
        return

    def new_cal_dipo_dislocations(self,
                                  in_tag="easy_Mid",
                                  input_s=0.0):
        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        self.set_lattce_constant(self._lattice_constant)
        self.set_element('W')

        n = 7
        m = 11
        atoms = self.set_bcc_convention([e1, e2, e3],
                                        (n, m, 1))  # z periodic 12
        movex = 0.0

        #  add tilt to the supercell #
        atoms = self.cut_half_atoms_new(atoms)
        supercell = atoms.get_cell()

        #  strain = \
        #  np.mat([[1.0, 0.0, -1./(3.)],
        #  [0.5, 1.0, (0.5 - 1./(6.))],
        #  [0.0, 0.0, 1.0]]);

        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.5, 1.0, 0.5],
                         [0.0, 0.0, 1.0]])

        supercell = strain * supercell
        print "new supercell is", supercell

        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])

        atoms_perf = atoms.copy()
        atoms = self.intro_dipole_screw_atoms(atoms, self._lattice_constant,
                                              movex, in_tag, input_s)

        #  atoms = self.intro_split_core(atoms);

        self.write_poscar(atoms_perf)
        os.system("cp POSCAR new_perf.vasp")
        self.write_poscar(atoms)
        os.system("cp POSCAR POSCAR_nofix")
        self.write_config_with_fix(atoms)
        self.prepare_dislocation_vasp_infiles()
        return

    def cal_dipo_dislocations(self):
        e1 = np.array([1., 1., -2.])
        e2 = np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        # 3, 5;   5, 9;  7, 11
        r, s = 7, 11
        v1 = r * e1
        v2 = 0.5 * (r * e1 + s * e2) + 0.5 * e3
        v3 = e3

        self.set_lattce_constant(self._lattice_constant)
        self.set_element('W')

        atoms = self.set_bcc_convention([v1, v2, v3], (r, 1, 1))
        print atoms.get_cell()
        atoms_perf = atoms.copy()
        atoms_perf = self.cut_half_atoms(atoms_perf)

        self.write_poscar(atoms_perf)
        os.system("cp POSCAR POSCAR_perf.vasp")

        atoms = self.intro_dipole_screw_atoms(atoms)
        atoms = self.cut_half_atoms(atoms)

        #  self.intro_dipole_screw()

        self.write_poscar(atoms)
        self.prepare_dislocation_vasp_infiles()
        return

    def prepare_dislocation_vasp_infiles(self):
        self.set_incar_type("dft2")
        self.write_incar()
        self.set_diff_kpoints([1, 2, 16])
        self.write_kpoints()

        self.set_nnodes(4)
        self.set_ppn(12)
        self.set_pbs_type('va')
        self.set_wall_time(200)
        self.set_job_title("w_easy_hard")
        self.set_main_job("mpirun vasp")

        self.write_pbs()
        return

    def cal_dislocations(self):
        m_direction = [[1, -2, 1], [1, 0, -1], [1, 1, 1]]
        # generate config #
        self.set_lattce_constant(self._lattice_constant)
        self.set_element('W')
        self.write_bcc_convention(m_direction, (5, 6, 3))

        #       Intro single screw          #
        self.intro_single_screw()
        return

    ############################################
    # generate config #
    ############################################

    def parepare_screw_dipo_with_solutes(self):
        sol_list = []
        for i in range(5):
            sol_list.append('%dc' % (i))
            sol_list.append('%dn' % (i))

        for i in range(len(sol_list)):

            dir_name = 'dir-%s' % (sol_list[i])
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)
            file_name = '%s_poscar.vasp' % (sol_list[i])
            with open(file_name, 'r') as fid:
                raw = fid.readlines()
                fid.close()
                print raw[5]
            raw[0] = "W   Ta\n"
            raw[5] = "W   Ta\n"

            #  with open("new", 'w') as fid:
            #  fid.writelines(raw);
            #  fid.close();

            #  os.system("mv  new  %s/POSCAR"%(dir_name));
            #  shutil.copy("POTCAR", dir_name);
            #  shutil.copy("KPOINTS", dir_name);
            #  shutil.copy("INCAR", dir_name);

            self.set_nnodes(4)
            self.set_pbs_type('va')
            self.set_wall_time(350)
            self.set_job_title("Ta_%s" % (dir_name))
            self.set_main_job("mpirun vasp")
            self.write_pbs()

            shutil.copy('va.pbs', dir_name)
            os.chdir(dir_name)
            os.system("qsub va.pbs")
            os.chdir(self.root_dir)
        return

    def test(self):
        in_tag = 'hard_split'
        s = 1
        self.new_cal_dipo_dislocations(in_tag, s)
        draw_dd = dd_map.dd_map_dislocation()
        draw_dd.draw_dd_map('init', in_tag, s)
        return

    def vasp_reaction_coordinate(self):
        in_tag = 'easy_hard'
        interp = 4
        delta = 1. / float(interp)
        draw_dd = dd_map.dd_map_dislocation()
        for i in range(interp + 1):
            s = i * delta
            dirname = 'dir-%d' % (i)

            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            os.chdir(dirname)

            self.new_cal_dipo_dislocations(in_tag, s)
            draw_dd.draw_dd_map('init', in_tag, s)
            os.system("cp init.png ../fig%d.png" % (i))
            os.system("cp ../POTCAR .")

            os.chdir(self.root_dir)
        return

    def loop_sub_jobs(self):
        for i in range(20):
            dirname = 'dir-{:03d}'.format(i)
            print dirname
            os.chdir(dirname)
            os.system("qsub va.pbs")
            os.chdir(os.pardir)
        return

    def collect_energy(self):
        flist = [0, 1, 2, 3, 5]
        engy, vol = np.zeros(len(flist)), np.zeros(len(flist))
        cnt = 0
        for i in flist:
            dirname = "dir-%d" % (i)
            os.chdir(dirname)
            (engy[cnt], vol[cnt]) = self.vasp_energy_stress_vol_quick()
            cnt += 1
            os.chdir(self.root_dir)
        print engy
        print vol
        return


usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--mtype", action="store",
                  type="string", dest="mtype", help="",
                  default="curv")
parser.add_option("-f", "--mfile", action="store",
                  type="string", dest="mfile",
                  default="./dummy.config.pair")

(options, args) = parser.parse_args()
OriDir = os.getcwd()

if __name__ == "__main__":
    N = vasp_dislocation()
    #  N.Intro_Screw_dipole()

    if options.mtype.lower() == 'new':
        N.new_cal_dipo_dislocations()

    if options.mtype.lower() == 'get':
        N.collect_energy()

    if options.mtype.lower() == 'path':
        N.vasp_reaction_coordinate()

    if options.mtype.lower() == 'sub':
        N.loop_sub_jobs()
