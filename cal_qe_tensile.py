#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-14 22:08:10
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-28 00:43:30
# encoding: utf-8

import os
import re
import shutil
import numpy as np
import gn_qe_inputs
import md_pot_data
import copy
from optparse import OptionParser

__version__ = 0.01
__author__ = 'Chaoming Yang'


class cal_ideal_tensile(object):

    def __init__(self,
                 Element,
                 AtomMass,
                 Percentage,
                 Lattice_constant,
                 Orientation,
                 Structure,
                 Kpoints,
                 Potential,
                 Num):

        self.pot = md_pot_data.qe_pot.vca_W75Re25

        self._tag = 'Opath'
        self._Num = Num
        if type(Element).__name__ == 'str':
            self._Element = Element
        elif type(Element).__name__ == 'list':
            self._Element1 = Element[0]
            self._Element2 = Element[1]
            self._Element = "%s%s" % (Element[0], Element[1])
            self._AtomMass1 = float(AtomMass[0])
            self._AtomMass2 = float(AtomMass[1])
            self._EffectiveMass = self._AtomMass1 * (1 - Percentage) + \
                self._AtomMass2 * Percentage

        if os.path.isdir("output"):
            print "mkdir output"
        else:
            os.mkdir("output")

        self._watchfile = "output/moniter_%d.txt" % (self._Num)
        self._outfile = "output/DATA"

        self._infile = "%s_%s.in" % (self._Element,
                                     self._Num)
        self._qeoutfile = "%s_%s.out" % (self._Element,
                                         self._Num)
        self._infile_backup = "%s_%s.in.back" % (self._Element,
                                                 self._Num)

        self._Lattice_constant = Lattice_constant
        self._Orientation = Orientation
        self._M = self.getM()
        self._Sij = self.getSij()
        self._Structure = Structure
        self._Kpoints = Kpoints
        self._addedstrain = np.mat([[0, 0, 0],
                                    [0, 0, 0],
                                    [0, 0, 0]],
                                   dtype="float")
        self._Stress_original = np.array([0, 0, 0, 0, 0, 0],
                                         dtype="float")
        self._Potential = Potential
        self._executable = "mpirun pw.x < %s > %s"\
            % (self._infile,
               self._qeoutfile)
        return

    def getM(self):
        orientation = self._Orientation
        x = np.array([0, 0, 0],
                     "float")
        y = np.array([0, 0, 0],
                     "float")
        for i in range(0, 3):
            x[i] = orientation[i]
            y[i] = orientation[3 + i]
        z = x.copy()
        z[0] = x[1] * y[2] - y[1] * x[2]
        z[1] = x[2] * y[0] - y[2] * x[0]
        z[2] = x[0] * y[1] - y[0] * x[1]
        xx, yy, zz = 0., 0., 0.
        for i in range(3):
            xx += x[i] * x[i]
            yy += y[i] * y[i]
            zz += z[i] * z[i]
        x /= np.sqrt(xx)
        y /= np.sqrt(yy)
        z /= np.sqrt(zz)
        M = np.mat([x, y, z], "float")
        return M

    def getSij(self):
        Element = self._Element1
        real_dn = r'-?(\d+\.\d*)'
        Element = '%s\s+' % (Element)
        print Element
        Object = Element + real_dn + \
            '\s+' + real_dn + '\s+' + real_dn + \
            '\s+' + real_dn + '\s+' + real_dn + '\s+' + real_dn
        print Object
        fid = open("cij_sij.txt", 'r')
        for line in fid:
            match = re.search(Object, line)
            if match:
                s11 = float(line.split()[4])
                s12 = float(line.split()[5])
                s44 = float(line.split()[6])
        fid.close()
        Sij = np.mat([[s11, s12, s12, 0, 0, 0],
                      [s12, s11, s12, 0, 0, 0],
                      [s12, s12, s11, 0, 0, 0],
                      [0, 0, 0, s44, 0, 0],
                      [0, 0, 0, 0, s44, 0],
                      [0, 0, 0, 0, 0, s44]], "float")
        return Sij

    def updateStrain(self,
                     InputStrain,
                     Correct_strain):
        addedstrain = np.zeros([3, 3],
                               "float")
        addedstrain[1, 1] = InputStrain[1]
        addedstrain[2, 2] = InputStrain[2]
        Correct_strain += addedstrain
        Correct_strain = np.mat(Correct_strain)
        return Correct_strain

    def get_data(self):
        with open(self._qeoutfile) as fid:
            Raw = fid.read()
            fid.close()
            print Raw
#            for line in Raw :
#                match = re.search(r"!\s*total energy\s*=\s*(-?\d*\.\d*)",line)
#                if match :
#                    print line
        findE = re.compile(r"!\s*total energy\s*=\s*(-?\d*\.\d*)")
        a = r"\s*total\s*stress\s*\(Ry\/bohr\*\*3\)\s*\(kbar\)\s*P=\s*(-?\d*\.\d*)\s*"
        b = r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*"
        c = r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*"
        d = r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)"
        M = a + b + c + d
        findStress = re.compile(M, re.DOTALL)
        Stress = findStress.findall(Raw)
        Energies = findE.findall(Raw)
        StressMatrix = np.zeros([3, 3],
                                dtype="float")
        for i in range(3):
            for j in range(3):
                StressMatrix[i, j] = 0.1 * \
                    float(Stress[-1][i * 6 + j + 4])  # Gpa
        print StressMatrix
        Stress6by1 = np.zeros([6, 1], "float")
        Stress6by1[0, 0] = StressMatrix[0, 0]
        Stress6by1[1, 0] = StressMatrix[1, 1]
        Stress6by1[2, 0] = StressMatrix[2, 2]
        Stress6by1[3, 0] = StressMatrix[1, 2]
        Stress6by1[4, 0] = StressMatrix[0, 2]
        Stress6by1[5, 0] = StressMatrix[0, 1]
        Stress6by1 = np.mat(Stress6by1)
        return (Energies[-1],
                Stress6by1)

    def output(self,
               Stress_original,
               delta,
               Correct_strain,
               Stress,
               Energy):

        with open(self._outfile, 'a') as ff:
            ff.write("delta %6.4f\t Sxx %10.8f\n"
                     % (delta, Stress_original[0]))
            print >> ff, "Stress original"
            print >> ff, Stress_original
            print >> ff, "calculated Stress"
            print >> ff, Stress
            print >> ff, "Energy is"
            print >> ff, Energy
            ff.close()
        if os.path.isfile(self._outfile):
            shutil.rmtree("./results")
        os.system("mv %s %s.in.%4.3f"
                  % (self._infile_backup, self._Element, delta))
        os.system("mv %s %s.out.%4.3f"
                  % (self._qeoutfile, self._Element, delta))
        return

    def gnInfile_T(self,
                   delta,
                   Correct_strain):
        original_strain = np.matrix([[1 + delta, 0, 0],
                                     [0, 1 - 0.1 * delta, 0],
                                     [0, 0, 1 + 0.1 * delta]],
                                    "float") + Correct_strain
        Base_vector = np.mat([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]], "float")
        Transformed_strain = np.mat(original_strain)
        Transposed_Base = Transformed_strain * Base_vector
        print np.mat(Transposed_Base)
        print self._Lattice_constant

        Transposed_Base = Transposed_Base * self._Lattice_constant
        AtomicPositions = np.mat([[0.0, 0.0, 0.0],
                                  [0.5, 0.5, 0.5]], "float")
        AtomicPositions = AtomicPositions * Transposed_Base

        return Transformed_strain

    def gnInfile_O(self,
                   delta,
                   Correct_strain):
        original_strain = np.matrix([[1 + delta, 0, 0],
                                     [0, 1 - 0.15 * delta, 0],
                                     [0, 0, 1 + 0.15 * delta]],
                                    "float") + Correct_strain

        Base_vector = np.mat([[1.0,     0,         0],
                              [0,    np.sqrt(2),   0],
                              [0,       0, np.sqrt(2)]], "float")

        Transformed_strain = np.mat(original_strain)

        Transposed_Base = Transformed_strain * Base_vector
        Transposed_Base = Transposed_Base * self._Lattice_constant

        AtomicPositions = np.mat([[0.0, 0.0, 0.0],
                                  [0.5, 0.0, 0.5],
                                  [0.5, 0.5, 0.0],
                                  [0.0, 0.5, 0.5]], "float")
        AtomicPositions = AtomicPositions * Transposed_Base

        with open(self._infile, 'w') as fid:
            fid.write("""
&control
calculation='relax',
prefix='%s',
tstress = .true.,
tprnfor = .true.,
outdir='./results',
pseudo_dir = './',
etot_conv_thr=1.0D-5,
forc_conv_thr=1.0D-4,
/
&system
    ibrav = 0, nat= %d, ntyp= 1,
    occupations='smearing',
    smearing='m-p',
    degauss=0.02D0,
    ecutwfc =45.0,
/
&electrons
conv_thr    = 1.D-6,
/
&ions
ion_dynamics='bfgs',
/
&cell
cell_dynamics='bfgs',
press_conv_thr=0.1D-0,
/
CELL_PARAMETERS {bohr}
%f  %f  %f
%f  %f  %f
%f  %f  %f
ATOMIC_SPECIES
%s  %f  %s
ATOMIC_POSITIONS {bohr}
%s  %f  %f  %f
%s  %f  %f  %f
%s  %f  %f  %f
%s  %f  %f  %f
K_POINTS automatic
%d %d %d  0  0  0
    """ % (self._Element,
                len(AtomicPositions),

                Transposed_Base[0, 0], Transposed_Base[
                    0, 1], Transposed_Base[0, 2],
                Transposed_Base[1, 0], Transposed_Base[
                    1, 1], Transposed_Base[1, 2],
                Transposed_Base[2, 0], Transposed_Base[
                    2, 1], Transposed_Base[2, 2],

                self._Element1,
                self._EffectiveMass,
                self._Potential,

                self._Element1,
                AtomicPositions[0, 0],
                AtomicPositions[0, 1],
                AtomicPositions[0, 2],

                self._Element1,
                AtomicPositions[1, 0],
                AtomicPositions[1, 1],
                AtomicPositions[1, 2],

                self._Element1,
                AtomicPositions[2, 0],
                AtomicPositions[2, 1],
                AtomicPositions[2, 2],

                self._Element1,
                AtomicPositions[3, 0],
                AtomicPositions[3, 1],
                AtomicPositions[3, 2],

                self._Kpoints[0],
                self._Kpoints[1],
                self._Kpoints[2],
           ))
        fid.close()
        return Transformed_strain

    def apply_modify_strain(self,
                            Stress,
                            StressPrime):
        StressPrime = float(StressPrime)
        Coeff =   0.18 + 0.02 * StressPrime ** 0.25 + \
            0.03 * StressPrime ** 0.19 + \
            0.05 * StressPrime ** 0.1
        if StressPrime > 20:
            Coeff =   0.20 + 0.01 * StressPrime ** 0.19 + \
                0.09 * StressPrime ** 0.1
        Strain = self._Sij * Stress * Coeff
        return Strain

    def loopcalculate(self):
        for i in range(6, 14):
            Correct_strain = np.mat([[0, 0, 0],
                                     [0, 0, 0],
                                     [0, 0, 0]], "float")

            delta = i * 0.02
            count = 0
            while True:
                if self._tag == 'Tpath':
                    if os.path.isfile(self._infile):
                        print "restart calculation"
                    else:
                        Transformed_strain = self.gnInfile_T(delta,
                                                             Correct_strain)
                elif self._tag == 'Opath':
                    if os.path.isfile(self._infile):
                        print "restart calculation"
                    else:
                        Transformed_strain = self.gnInfile_O(delta,
                                                             Correct_strain)
                else:
                    print "please specify Opath or Tpath"

                os.system("%s" % (self._executable))
                # back up infile
                shutil.move(self._infile, self._infile_backup)

                Energy, Stress_original = self.get_data()

                Stress = copy.deepcopy(Stress_original)
                Stress[0, 0] = 0.0
                Stress_abs = np.abs(Stress)
                StressPrime = np.max([Stress_abs[1, 0],
                                      Stress_abs[2, 0]])

                if StressPrime < 1.0:
                    self.output(Stress_original,
                                delta,
                                Correct_strain,
                                Stress,
                                Energy)
                    break
                else:
                    Strain = self.apply_modify_strain(Stress,
                                                      StressPrime)

                    Correct_strain = self.updateStrain(Strain,
                                                       Correct_strain)
                    with open(self._watchfile, 'a') as fout:
                        print >> fout, "loop time",  count
                        print >> fout, "StressPrime", StressPrime
                        print >> fout, "original Stress", Stress_original
                        print >> fout, "Correct_strain", Correct_strain
                        count += 1
                        fout.close()
        return

    def LoopCalMoreAccurate(self):
        RootDir = os.getcwd()
        for i in range(5, 13):
            count = 0
            Correct_strain = np.mat([[0, 0, 0],
                                     [0, 0, 0],
                                     [0, 0, 0]], "float")
            delta = i * 0.02
            Dir = "Dir-%4.3f" % (delta)
            originInfile = "WRe.in.%4.3f" % (delta)
            print originInfile

            os.mkdir(Dir)
            if os.path.isfile(originInfile):
                shutil.copy(originInfile, Dir)
            shutil.copy(self._Potential, Dir)
            shutil.copy("cij_sij.txt", Dir)
            os.chdir(Dir)
            os.mkdir("output")

            while True:
                if self._tag == 'Tpath':
                    if os.path.isfile(originInfile):
                        print "restart calculation"
                        shutil.move(originInfile, self._infile)
                    else:
                        Transformed_strain = self.gnInfile_T(delta,
                                                             Correct_strain)
                elif self._tag == 'Opath':
                    if os.path.isfile(originInfile):
                        print "restart calculation"
                        shutil.move(originInfile, self._infile)
                    else:
                        Transformed_strain = self.gnInfile_O(delta,
                                                             Correct_strain)
                else:
                    print "please specify Opath or Tpath"

                os.system("%s" % (self._executable))
                # back up infile
                shutil.move(self._infile, self._infile_backup)

                Energy, Stress_original = self.get_data()

                Stress = copy.deepcopy(Stress_original)
                Stress[0, 0] = 0.0
                Stress_abs = np.abs(Stress)
                StressPrime = np.max([Stress_abs[1, 0],
                                      Stress_abs[2, 0]])

                if StressPrime < 0.05:
                    self.output(Stress_original,
                                delta,
                                Correct_strain,
                                Stress,
                                Energy)
                    break
                else:
                    Strain = self.apply_modify_strain(Stress,
                                                      StressPrime)

                    Correct_strain = self.updateStrain(Strain,
                                                       Correct_strain)
                    with open(self._watchfile, 'a') as fout:
                        print >> fout, "loop time",  count
                        print >> fout, "StressPrime", StressPrime
                        print >> fout, "original Stress", Stress_original
                        print >> fout, "Correct_strain", Correct_strain
                        count += 1
                        fout.close()
            os.chdir(RootDir)
        return


def loop_potentials():
    RootDir = os.getcwd()
    Potentialist = ['WRe.0-05.fhi.UPF',
                    'WRe.0-1.fhi.UPF',
                    'WRe.0-15.fhi.UPF',
                    'WRe.0-2.fhi.UPF',
                    'WRe.0-25.fhi.UPF']

    LatticeList = [5.98753,
                   5.98281,
                   5.97809,
                   5.97337,
                   5.92]

    for i in range(1, 2):
        Potential = Potentialist[i]
        Lattice_constant = LatticeList[i]

        DirName = "Dir-%s" % (Potential)
        os.mkdir(DirName)
        shutil.copy(Potential, DirName)
        shutil.copy('cij_sij.txt', DirName)
        os.chdir(DirName)
        os.mkdir('results')

        M = QEtensile(Element=['W', 'Re'],
                      AtomMass=[183.84, 186.207],
                      Percentage=(i + 1) * 0.05,
                      Lattice_constant=Lattice_constant,
                      Orientation=[1, 0, 0, 0, 1, 0],
                      Structure='bcc',
                      Kpoints=[25, 25, 25],
                      Potential=Potential)
        M.loopcalculate()

        FileOut = "DATA-%s" % (Potential)
        FileDest = os.path.join(RootDir, FileOut)
        shutil.copy("DATA", FileDest)

        os.chdir(RootDir)
    return


def calculate(PotentialN):
    Potentialist = ['WRe.0-05.fhi.UPF',
                    'WRe.0-10.fhi.UPF',
                    'WRe.0-15.fhi.UPF',
                    'WRe.0-20.fhi.UPF',
                    'WRe.0-25.fhi.UPF',
                    'WRe.0-50.fhi.UPF']

    LatticeList = [5.98753,
                   5.98281,
                   5.97809,
                   5.97337,
                   5.96924,
                   5.94977]

    M = QEtensile(Element=['W', 'Re'],
                  AtomMass=[183.84, 186.207],
                  Percentage=(PotentialN + 1) * 0.05,
                  Lattice_constant=LatticeList[PotentialN],
                  Orientation=[1, 0, 0, 0, 1, 0],
                  Structure='bcc',
                  Kpoints=[25, 20, 20],
                  Potential=Potentialist[PotentialN],
                  Num=0)
    M.LoopCalMoreAccurate()
    return

usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--structure",
                  dest="structure",
                  help="fcc bcc or hcp")

parser.add_option("-d", "--direction",
                  dest="orientation",
                  nargs=6, type="int",
                  help=" x y z yz xz xy",
                  default="1")

parser.add_option("-e",
                  "--element",
                  dest="Element",
                  help="type of Element",
                  default="Mo")

parser.add_option("-t",
                  "--mtype",
                  action="store",
                  type="string",
                  dest="mtype",
                  help="",
                  default="prp_r")

parser.add_option("-l", "--lattice_constant",
                  dest="lattice_constant",
                  type="float",
                  help="Lattice constant",
                  default=3.15)
(options, args) = parser.parse_args()

if __name__ == "__main__":
    calculate(4)
    if options.mtype.lower() == 'prepqe':
        calculation()

    if options.mtype.lower() == 'clc':
        loop_clc_output()
