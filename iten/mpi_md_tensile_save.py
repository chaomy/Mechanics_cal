#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-05 23:59:58


import copy
import os
import numpy as np
from multiprocessing import Pool
import shutil
import sys

try:
    import get_data
    import output_data
    import gn_lmp_infile

except ImportError:
    print("error during import")

__version__ = 0.01
__author__ = 'Chaoming Yang'
__all__ = ['md_tensile']


def unwrap_self_run_lammps(arg, **kwarg):
    return md_loop_tensile.lammps_job(*arg, **kwarg)


class md_tensile(get_data.get_data,
                 output_data.output_data,
                 gn_lmp_infile.gn_md_infile):

    def __init__(self, element,
                 lattice_constant, size,
                 orientation, structure):

        get_data.get_data.__init__(self)
        output_data.output_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)

        self.element = element
        self.lattice_constant = lattice_constant
        self.set_lattice_constant()

        self.size = size
        self.orientation = orientation
        self.M = self.get_M()
        self.Sij = self.get_Sij()
        self.structure = structure

        self.lx, self.ly, self.lz = 0, 0, 0
        self.addedstrain = np.mat([[0, 0, 0],
                                   [0, 0, 0],
                                   [0, 0, 0]], "float")

        self.stress_original = np.array([0, 0, 0, 0, 0, 0],
                                        dtype="float")

    def update_strain(self, Inputstrain,
                      stress, Correct_strain):
        addedstrain = copy.deepcopy(self.addedstrain)

        if abs(stress[1]) > self._stress_ThrValue:
            addedstrain[1, 1] = Inputstrain[1]
            Correct_strain[1, 1] += addedstrain[1, 1]
        if abs(stress[2]) > self._stress_ThrValue:
            addedstrain[2, 2] = Inputstrain[2]
            Correct_strain[2, 2] += addedstrain[2, 2]
        return Correct_strain

    def gn_bcc_tpath(self,
                     delta,
                     Correct_strain):

        element = self.element
        M = self.get_M()
        lattice_constant = self.lattice_constant
        size = self.size
        M = self.M

        #  original_strain = np.matrix([[1 + delta, 0,    0],
        #  [0, 1 + 0.2 * delta, 0],
        #  [0, 0, 1 - 0.2 * delta]],
        #  "float") + Correct_strain

        Transformed_strain = M.transpose() * self.strainmtx * M
        Base_vector = np.matrix([[1, 0, 0],
                                 [0, 1, 0],
                                 [0, 0, 1]])
        Transposed_Base = Transformed_strain * Base_vector

        xsize = size[0]
        ysize = size[1]
        zsize = size[2]
        a = 1. / xsize
        b = 1. / ysize
        c = 1. / zsize

        Lengthx = Transposed_Base[0, :] * Transposed_Base[0, :].transpose()
        Lengthy = Transposed_Base[1, :] * Transposed_Base[1, :].transpose()
        Lengthz = Transposed_Base[2, :] * Transposed_Base[2, :].transpose()
        Lengthx = np.sqrt(Lengthx)
        Lengthx *= lattice_constant
        Lengthy = np.sqrt(Lengthy)
        Lengthy *= lattice_constant
        Lengthz = np.sqrt(Lengthz)
        Lengthz *= lattice_constant

        xlo, xhi = 0.0, Lengthx * xsize
        ylo, yhi = 0.0, Lengthy * ysize
        zlo, zhi = 0.0, Lengthz * zsize

        AtomNumber = int(xsize * ysize * zsize) * 2
        XDirection, YDirection, ZDirection = [], [], []

        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XDirection.append((x) * a * Lengthx)
                    YDirection.append((y) * b * Lengthy)
                    ZDirection.append((z) * c * Lengthz)

                    XDirection.append((x + 0.5) * a * Lengthx)
                    YDirection.append((y + 0.5) * b * Lengthy)
                    ZDirection.append((z + 0.5) * c * Lengthz)

        filename = "lattice.txt"  # %(element,lattice_constant)

        with open(filename, mode="w") as fout:
            fout.write("#bcc_structure_input for lammps\n")
            fout.write("\n")
            fout.write("%d atoms\n" % (AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (xlo, xhi))
            fout.write("%f\t%f ylo yhi\n" % (ylo, yhi))
            fout.write("%f\t%f zlo zhi\n" % (zlo, zhi))
            fout.write("%8.5f %8.5f %8.5f xy xz yz" % (0, 0, 0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                           % (i + 1,
                              XDirection[i],
                              YDirection[i],
                              ZDirection[i]))
        fout.close()
        os.system("cp lattice.txt lattice_%5.3f.txt" % (delta))
        ########### generate cfg #########################
        XXDirection, YYDirection, ZZDirection = [], [], []
        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XXDirection.append((x) * a)
                    YYDirection.append((y) * b)
                    ZZDirection.append((z) * c)

                    XXDirection.append((x + 0.5) * a)
                    YYDirection.append((y + 0.5) * b)
                    ZZDirection.append((z + 0.5) * c)

        curdir = os.getcwd()
        os.chdir("mycfg")
        filename = "cfg_%5.4f.cfg" % (delta)
        with open(filename, mode="w") as fout:
            fout.write("Number of particles = %d\n" % (AtomNumber))
            fout.write("""A = 1 Angstrom
H0(1,1) = %f A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = %f A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = %f A
Transform(1,1) = 1
Transform(1,2) = %f
Transform(1,3) = %f
Transform(2,1) = %f
Transform(2,2) = 1
Transform(2,3) = %f
Transform(3,1) = %f
Transform(3,2) = %f
Transform(3,3) = 1
    """ % (Transposed_Base[0, 0] * xsize * Lengthx,
                Transposed_Base[1, 1] * ysize * Lengthy,
                Transposed_Base[2, 2] * zsize * Lengthz,
                Transposed_Base[0, 1] * xsize,
                Transposed_Base[0, 2] * xsize,
                Transposed_Base[1, 0] * xsize,
                Transposed_Base[1, 2] * xsize,
                Transposed_Base[2, 0] * xsize,
                Transposed_Base[2, 1] * xsize))

            for i in range(AtomNumber):
                fout.write("95.94\t%s\t%7.7f\t %7.7f\t %7.7f\t0\t0\t0\n"
                           % (element,
                              XXDirection[i],
                              YYDirection[i],
                              ZZDirection[i]))
        fout.close()
        os.chdir(curdir)
        return

    def gn_bcc_opath(self,
                     delta,
                     Correct_strain):
        element = self.element
        M = self.get_M()
        lattice_constant = self.lattice_constant
        size = self.size
        M = self.M

        #  original_strain = np.matrix([[1 + delta, 0, 0],
        #  [0, 1 + 0.2 * delta, 0],
        #  [0, 0, 1 - 0.2 * delta]],
        #  "float") + Correct_strain

        Transformed_strain = M.transpose() * self.strainmtx * M
        Base_vector = np.matrix([[1, 0, 0],
                                 [0, np.sqrt(2), 0],
                                 [0, 0, np.sqrt(2)]])

        Transposed_Base = Transformed_strain * Base_vector
        xsize = size[0]
        ysize = size[1]
        zsize = size[2]
        a = 1. / xsize
        b = 1. / ysize
        c = 1. / zsize
        Lengthx = Transposed_Base[0, :] * Transposed_Base[0, :].transpose()
        Lengthy = Transposed_Base[1, :] * Transposed_Base[1, :].transpose()
        Lengthz = Transposed_Base[2, :] * Transposed_Base[2, :].transpose()

        Lengthx = np.sqrt(Lengthx)
        Lengthx *= lattice_constant
        Lengthy = np.sqrt(Lengthy)
        Lengthy *= lattice_constant
        Lengthz = np.sqrt(Lengthz)
        Lengthz *= lattice_constant

        xlo, xhi = 0.0, Lengthx * xsize
        ylo, yhi = 0.0, Lengthy * ysize
        zlo, zhi = 0.0, Lengthz * zsize

        AtomNumber = int(xsize * ysize * zsize) * 4
        XDirection, YDirection, ZDirection = [], [], []

        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XDirection.append((x) * a * Lengthx)
                    YDirection.append((y) * b * Lengthy)
                    ZDirection.append((z) * c * Lengthz)

                    XDirection.append((x) * a * Lengthx)
                    YDirection.append((y + 0.5) * b * Lengthy)
                    ZDirection.append((z + 0.5) * c * Lengthz)

                    XDirection.append((x + 0.5) * a * Lengthx)
                    YDirection.append((y) * b * Lengthy)
                    ZDirection.append((z + 0.5) * c * Lengthz)

                    XDirection.append((x + 0.5) * a * Lengthx)
                    YDirection.append((y + 0.5) * b * Lengthy)
                    ZDirection.append((z) * c * Lengthz)

        filename = "lattice.txt"  # %(element,lattice_constant)
        with open(filename, mode="w") as fout:
            fout.write("#bcc_structure_input for lammps\n")
            fout.write("\n")
            fout.write("%d atoms\n" % (AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (xlo, xhi))
            fout.write("%f\t%f ylo yhi\n" % (ylo, yhi))
            fout.write("%f\t%f zlo zhi\n" % (zlo, zhi))
            fout.write("%8.5f %8.5f %8.5f xy xz yz" % (0, 0, 0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                           % (i + 1,
                              XDirection[i],
                              YDirection[i],
                              ZDirection[i]))
        fout.close()
        os.system("cp lattice.txt lattice_%5.3f.txt" % (delta))
        # generate cfg #
        XXDirection, YYDirection, ZZDirection = [], [], []
        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XXDirection.append((x) * a)
                    YYDirection.append((y) * b)
                    ZZDirection.append((z) * c)

                    XXDirection.append((x) * a)
                    YYDirection.append((y + 0.5) * b)
                    ZZDirection.append((z + 0.5) * c)

                    XXDirection.append((x + 0.5) * a)
                    YYDirection.append((y) * b)
                    ZZDirection.append((z + 0.5) * c)

                    XXDirection.append((x + 0.5) * a)
                    YYDirection.append((y + 0.5) * b)
                    ZZDirection.append((z) * c)
        curdir = os.getcwd()
        os.chdir("mycfg")
        filename = "cfg_%5.4f.cfg" % (delta)
        with open(filename, mode="w") as fout:
            fout.write("Number of particles = %d\n" % (AtomNumber))
            fout.write("""A = 1 Angstrom
H0(1,1) = %f A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = %f A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = %f A
Transform(1,1) = 1
Transform(1,2) = %f
Transform(1,3) = %f
Transform(2,1) = %f
Transform(2,2) = 1
Transform(2,3) = %f
Transform(3,1) = %f
Transform(3,2) = %f
Transform(3,3) = 1
""" % (Transposed_Base[0, 0] * xsize * Lengthx,
                Transposed_Base[1, 1] * ysize * Lengthy,
                Transposed_Base[2, 2] * zsize * Lengthz,
                Transposed_Base[0, 1] * xsize,
                Transposed_Base[0, 2] * xsize,
                Transposed_Base[1, 0] * xsize,
                Transposed_Base[1, 2] * xsize,
                Transposed_Base[2, 0] * xsize,
                Transposed_Base[2, 1] * xsize))
            for i in range(AtomNumber):
                fout.write("95.94\t%s\t%7.7f\t %7.7f\t %7.7f\t0\t0\t0\n"
                           % (element,
                              XXDirection[i],
                              YYDirection[i],
                              ZZDirection[i]))
        fout.close()
        os.chdir(curdir)
        return


class md_loop_tensile(md_tensile):

    def __init__(self,
                 element='Nb',
                 lattice_constant=3.30,
                 size=[1, 1, 1],
                 orientation=[[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]],
                 structure='bcc',
                 in_potential='dummy.lammps.ADP'):

        self._flux_exe = 'lmp_mpi < in.stat_tensile'

        self._looptime = 10
        self._increment = 0.02
        self._stress_ThrValue = 0.05

        self._tensile_potential = in_potential
        self._element = element
        self._lattice_constant = lattice_constant
        self._size = size
        self._orientation = orientation

        self._structure = structure
        self.root_dir = os.getcwd()
        return

    def set_orientation(self, orientation):
        self._orientation = orientation
        return

    def set_lattce_constant(self, lattice_constant):
        self._lattice_constant = lattice_constant
        return

    def set_loop_time(self, loop_time):
        self._looptime = loop_time
        return

    def set_potential(self, potential):
        self._tensile_potential = potential
        return

    def get_strain(self,
                   Sij,
                   stress,
                   stressPrime):
        tag = "constant"
        strain = np.zeros(6).transpose()
        if tag == "constant":
            coeff = 1.0
            strain[1] = Sij[0, 1] * stress[1] * coeff
            strain[2] = Sij[0, 2] * stress[2] * coeff
        else:
            stressPrime = float(stressPrime)
            coeff = 0.06 + 0.06 * stressPrime ** 0.20 + \
                0.03 * stressPrime ** 0.19 + \
                0.01 * stressPrime ** 0.1
            if stressPrime > 20:
                coeff = 0.08 + 0.01 * stressPrime ** 0.19 + \
                    0.09 * stressPrime ** 0.1
            strain = -Sij * stress * coeff
        return strain

    def lammps_job(self,
                   option,
                   job):

        if option == 'TP' or option == 'tpath':
            directory = "TP_%s_%d" % (self._element,
                                      job)

            output_data_file = os.path.join(self.root_dir, 'TP-DATA')
        elif option == 'OP' or option == 'opath':
            directory = "OP_%s_%d" % (self._element,
                                      job)

            output_data_file = os.path.join(self.root_dir, 'OP-DATA')

        if os.path.isdir(directory):
            shutil.rmtree(directory)

        os.mkdir(directory)
        shutil.copy("./in.stat_tensile", directory)
        shutil.copy(self._tensile_potential, directory)
        os.system("cp ~/src/Data_process/cij_sij.txt %s" % (directory)) 
        os.chdir(directory)
        self.gn_md_tensile(potential_file=self._tensile_potential,
                           element=self._element)
        os.mkdir("cfg")
        os.mkdir("mycfg")
        os.mkdir("dump")
        os.mkdir("restart")
        os.mkdir("xyz")

        if os.path.isfile(output_data_file):
            os.system(": > %s" % (output_data_file))
        # main function #
        Correct_strain = np.mat([[0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]], "float")
        md_tensile.__init__(self,
                            self._element,
                            self._lattice_constant,
                            self._size,
                            self._orientation,
                            self._structure)

        delta = job * 0.02
        count = 0
        coeff = 1.0

        self.strainmtx = np.matrix([[1 + delta, 0, 0],
                                    [0, 1 + 0.2 * delta, 0],
                                    [0, 0, 1 - 0.2 * delta]],
                                   "float") + Correct_strain

        strain = np.zeros(6).transpose()
        while True:
            if option == 'TP':
                self.gn_bcc_tpath(delta, Correct_strain)

            elif option == 'OP':
                self.gn_bcc_opath(delta, Correct_strain)

            #  os.system("cat lattice.txt >> backup.txt")
            os.system("%s > Log_MD" % (self._flux_exe))

            # lx, ly, lz #
            self.lx, self.ly, self.lz = self.md_get_lx_ly_lz_from_Log("Log_MD")

            self.stress_original = self.md_get_stress()
            stress = copy.deepcopy(self.stress_original)

            stress[0] = 0.0
            stress_abs = np.abs(stress)
            stressPrime = np.max(stress_abs)

            if stressPrime < self._stress_ThrValue:
                self.output_md_tensile(delta)
                break
            elif count > 500:
                self.output_md_tensile(delta)
                break
            else:
                # update strain #
                print("Sij = ", self.Sij[0, 1], self.Sij[0, 2])

                if abs(stress[1]) > self._stress_ThrValue:
                    strain[1] = self.Sij[0, 1] * stress[1] * coeff
                    self.strainmtx[1, 1] += strain[1]

                if abs(stress[2]) > self._stress_ThrValue:
                    strain[2] = self.Sij[0, 2] * stress[2] * coeff
                    self.strainmtx[2, 2] += strain[2]

                #  Correct_strain = self.update_strain(strain,
                    #  stress,
                    #  Correct_strain)

                with open('monitor.txt', 'a') as fid:
                    print("delta ", delta, file=fid)
                    print("Run times", count, file=fid)
                    print("stress_original", self.stress_original, file=fid)
                    print("stressPrime ", stressPrime, file=fid)
                    fid.close()

                count += 1
        os.system("cat DATA >> %s" % (output_data_file))
        os.chdir(self.root_dir)
        return (delta, self.stress_original[0, 0])

    def cal_md_tensile(self,
                       options):
        pool = Pool(processes=self._looptime)
        List = np.arange(self._looptime)
        results = pool.map(unwrap_self_run_lammps,
                           list(zip([self] * len(List),
                               [options] * len(List),
                               List)))
        strain = []
        stress = []
        for i in range(self._looptime):
            strain.append(results[i][0])
            stress.append(results[i][1])
        return (strain, stress)

    def output_log(self, strain, stresstp, stressop):
        with open("tensile.log", "w") as fid:
            for i in range(len(strain)):
                fid.write("%8.5f  %10.5f  %10.5f\n"
                          % (strain[i],
                             stresstp[i],
                             stressop[i]))
            fid.close()
        #  os.system("rm -rf TP_%s_*" % (self._element))
        #  os.system("rm -rf OP_%s_*" % (self._element))
        return

    def set_lattice_constant(self, opt="given"):
        if opt == "given":
            self.lattice_constant = float(sys.argv[1])
        return


if __name__ == "__main__":
    makePot = "dummy.lammps.ADP"
    givenPot = "Nb.eam.alloy.webarchive"
    M = md_loop_tensile(element='Nb',
                        lattice_constant=3.34349187552591,
                        size=[1, 1, 1],
                        orientation=[1, 0, 0, 0, 1, 0],
                        structure='bcc',
                        in_potential=makePot)
    (strain, stress) = M.cal_md_tensile('TP')
    # (strain2, stress2) = M.cal_md_tensile('OP')
    # M.output_log(strain, stress, stress2)
