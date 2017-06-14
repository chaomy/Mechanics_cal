#!/usr/bin/env python

###################################################################
#
# File Name : ./cal_vasp_cij.py
#
###################################################################
#
# Purpose : cal cij by vasp
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import glob
from optparse import OptionParser
import os
import numpy as np
from scipy.optimize import leastsq

try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import output_data

except ImportError:
    print "error during import"


class cal_cij(gn_config.bcc,
              gn_config.fcc,
              gn_config.hcp,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              output_data.output_data):

    def __init__(self,
                 in_structure='bcc',
                 lattice_constant=3.1711,
                 in_element='W',
                 in_kpoints=[31, 31, 31]):

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        output_data.output_data.__init__(self)

        self.cij_unit_delta = 0.005
        self.cij_loop_times = 10
        self.volume = None
        self.energy0 = None

        self._cij_lattice_constant = lattice_constant
        self._cij_element = in_element
        self.cij_type_list = ['c11', 'c12', 'c44']  # default
        self.cij_type = 'c11'
        self._surface_kpoints = in_kpoints
        self._structure = in_structure

        if self.structure == 'bcc':
            gn_config.bcc.__init__(self)
        elif self.structure == 'fcc':
            gn_config.fcc.__init__(self)
        elif self.structure == 'hcp':
            gn_config.hcp.__init__(self)
        return

    def fit_para(self, delta_list, energy_list):
        xdata = delta_list
        ydata = energy_list - self.energy0

        def residuals(p):
            a, c = p
            return ydata - (a * (0.5 * xdata**2) + c)
        r = leastsq(residuals, [1, 0])
        a, c = r[0]
        a = a / self.volume * self.ev_angstrom3_to_GPa
        return a

    def obtain_cij(self, opt='np'):

        self.set_volume_energy0()
        if opt is 'np':
            raw1 = np.loadtxt("data_c11.txt")
            raw2 = np.loadtxt("data_c12.txt")
            raw3 = np.loadtxt("data_c44.txt")
        else:
            delta_list, energy_list = self.get_delta_energy("data_c11.txt")

        delta_list, energy_list = raw1[:, 0], raw1[:, 1]
        c11_plus_c12 = self.fit_para(delta_list, energy_list)

        delta_list, energy_list = raw2[:, 0], raw2[:, 1]
        c11_minus_c12 = self.fit_para(delta_list, energy_list)

        delta_list, energy_list = raw3[:, 0], raw3[:, 1]
        c44 = self.fit_para(delta_list, energy_list)

        c11 = (0.111111111111111 * c11_plus_c12 +
               0.333333333333333 * c11_minus_c12)
        c12 = (0.111111111111111 * c11_plus_c12 -
               0.166666666666667 * c11_minus_c12)
        c44 = 0.25 * c44

        with open("cij.dat", 'w') as fout:
            fout.write("C11\t%f\t\nC12\t%f\t\nC44\t%f\t\n" % (c11, c12, c44))
            fout.close()
        return

    def obtain_cij_old(self):
        self.volume = 18.30
        self.energy0 = -10.09407662

        data = np.loadtxt("c11_summary_35.00")
        delta_list = np.array(data[:, 0])
        energy_list = np.array(data[:, 1])
        c11_plus_c12 = self.fit_para(delta_list, energy_list)

        data = np.loadtxt("c12_summary_35.00")
        delta_list = np.array(data[:, 0])
        energy_list = np.array(data[:, 1])
        c11_minus_c12 = self.fit_para(delta_list, energy_list)

        data = np.loadtxt("c44_summary_35.00")
        delta_list = np.array(data[:, 0])
        energy_list = np.array(data[:, 1])
        c44 = self.fit_para(delta_list, energy_list)

        c11 = (0.111111111111111 * c11_plus_c12 +
               0.333333333333333 * c11_minus_c12)
        c12 = (0.111111111111111 * c11_plus_c12 -
               0.166666666666667 * c11_minus_c12)
        c44 = 0.25 * c44

        print c11, c12, c44
        return

    def set_volume_energy0(self):
        raw = np.loadtxt("equilibrium.txt")
        self.volume = np.average(raw[:, 0])
        self.energy0 = np.average(raw[:, 1])
        return

    def set_cij_type(self, cij_type):
        self.cij_type = cij_type
        return

    def write_cij_poscar(self, delta):
        self.set_lattce_constant(self._cij_lattice_constant)
        self.write_bcc_primitive_with_strain(delta,
                                             self.cij_type,
                                             (1, 1, 1))
        return

    def loop_prepare_cij(self):
        for i in range(len(self.cij_type_list)):
            current_type = self.cij_type_list[i]
            self.set_cij_type(current_type)

            for j in range(-self.cij_loop_times, self.cij_loop_times):
                delta = self.cij_unit_delta * j
                if j >= 0:
                    dir_name = "dir-%s-delta-%03d" % (current_type, j)
                else:
                    dir_name = "dir-%s-delta-n%03d" % (current_type, np.abs(j))
                if not os.path.isdir(dir_name):
                    os.mkdir(dir_name)
                os.chdir(dir_name)

                self.write_cij_poscar(delta)
                self.prepare_vasp_inputs(dir_name)
                os.chdir(self.root_dir)
        return

    ###################################################################
    # continue the calculation due to the limit of walltime
    ###################################################################
    def cal_cij_continue(self, tag='modify'):
        dirlist = glob.glob("dir-c*-delta-*")
        for i in range(len(dirlist)):
            mydir = dirlist[i]
            os.chdir(mydir)
            os.system("cp CONTCAR POSCAR")

            if tag is 'cnt':
                # rewrite the va.pbs #
                self.set_pbs_type('va')
                self.set_wall_time(50)
                self.set_job_title(mydir[6:])
                self.set_nnodes(1)
                self.set_ppn(12)
                self.set_main_job("mpirun vasp")
                self.write_pbs(od=None)

                # sub the job
                os.system("qsub va.pbs")

            elif tag is 'modify':
                #  os.system("rm OUTCAR_cpy")
                #  os.system("head -n 500 OUTCAR >> OUTCAR_cpy")
                #  os.system("tail -n 1000 OUTCAR >> OUTCAR_cpy")
                #  os.system("rm OUTCAR")
                os.system("mv OUTCAR_cpy OUTCAR")
            os.chdir(self.root_dir)
        return

    def ouput_data_cij(self):
        if os.path.isfile("cij.dat"):
            os.system("mv cij.dat cij.dat")

        for i in range(len(self.cij_type_list)):
            current_type = self.cij_type_list[i]
            self.set_cij_type(current_type)
            out_file_name = "data_%s.txt" % (current_type)

            for j in range(-self.cij_loop_times, self.cij_loop_times):
                delta = self.cij_unit_delta * j
                if j >= 0:
                    dir_name = "dir-%s-delta-%03d" % (current_type, j)
                else:
                    dir_name = "dir-%s-delta-n%03d" % (current_type, np.abs(j))
                os.chdir(dir_name)
                print "i am in ", dir_name
                (energy, volume) = self.vasp_energy_stress_vol_quick()
                os.chdir(self.root_dir)

                self.output_delta_energy(delta, energy, file_name=out_file_name)
                if j == 0:
                    self.output_equilibrium(energy=energy,
                                            volume=volume)

        #  fout = open("%s_summary"%cname,"w")
            #  fout.write("%22.16f %22.16f\n"%(delta_increment*i, energy))
        #  answer = fit_para(cname,E0)
        return

    def prepare_vasp_inputs(self, dir_name):
        self.set_incar_type('dft')
        self.set_accuracy(1e-4)
        self.write_incar()

        self.set_diff_kpoints([31, 31, 31])
        self.set_intype('gamma')
        self.write_kpoints()

        self.set_pbs_type('va')
        self.set_wall_time(30)
        self.set_job_title(dir_name[6:])
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=None)

        os.system("cp ../../POTCAR .")
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

if __name__ == "__main__":

    Job = cal_cij(in_structure='bcc',
                  lattice_constant=3.322404,
                  in_element='Nb',
                  in_kpoints=[31, 31, 31])

    if options.mtype.lower() == 'prep':
        Job.loop_prepare_cij()

    if options.mtype.lower() == 'sub':
        Job.loop_sub_jobs()

    if options.mtype.lower() == 'clc':
        Job.obtain_cij()

    if options.mtype.lower() == 'cnt':
        Job.cal_cij_continue()

    if options.mtype.lower() == 'data':
        Job.ouput_data_cij()

    #  A.set_volume_energy0()
    #  A.obtain_cij_old()

    #  c12 = 0.25 * C11plus2C12 + 0.75 * C11minusC12
    #  c11 = 0.25 * C11plus2C12 - 0.25 * C11minusC12
    #  c44 = float(get_cxx("c44",Lattice_constant,delta_increment,nt,E0)) * (1./4.) * inv_V0  * unit_conversion1
