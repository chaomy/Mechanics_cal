#!/usr/bin/env python
# encoding: utf-8
import os
import numpy as np
from optparse import OptionParser
import copy
__version__ = 0.01
__author__ = 'Chaoming Yang'
__all__ = ['VAtension']


class VAtension(object):
    def __init__(self, Element, lattice_constant, orientation, structure):
        self.Element = Element
        self.lattice_constant = lattice_constant
        self.orientation = orientation
        self.M = self.getM()
        self.Sij = self.getSij()
        self.structure = structure
        self.addedstrain = np.mat([[0, 0, 0], [0, 0, 0], [0, 0, 0]], "float")
        self.__Stress_original = np.array([0, 0, 0, 0, 0, 0], dtype="float")
        return

    def getM(self):
        orientation = self.orientation
        x = np.array([0, 0, 0], "float")
        y = np.array([0, 0, 0], "float")
        for i in range(0, 3):
            x[i] = orientation[i]
            y[i] = orientation[3 + i]
        z = x.copy()
        z[0] = x[1] * y[2] - y[1] * x[2]
        z[1] = x[2] * y[0] - y[2] * x[0]
        z[2] = x[0] * y[1] - y[0] * x[1]

        ######## unify ##############
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
        Element = self.Element
        import re
        real_dn = r'-?(\d+\.\d*)'
        Element = '%s\s+' % (Element)
        print(Element)
        Object = Element + real_dn + '\s+' + real_dn + '\s+' + \
            real_dn + '\s+' + real_dn + '\s+' + real_dn + '\s+' + real_dn
        print(Object)
        fid = open("./cij_sij.txt", 'r')
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

    def grabStress(self):
        import re

        def get_data(expr):
            os.system("%s>eout" % (expr))
            fin = file("eout", 'r')
            line = fin.readline().split()
            return (float(line[0]))
        Energy = get_data("tail -n1 OSZICAR |awk '{print $5}'")
        Stress = np.arange(6, dtype="float")
        Stress.shape = ([6, 1])
        with open("OUTCAR", 'r') as fid:
            real_N = real_N = r"(-?\d*\.\d{5})"
            In = r"\s*"
            for line in fid:
                match = re.search(r"^\s*in kB\s*" + real_N + In + real_N + In +
                                  real_N + In + real_N + In + real_N + In + real_N, line)
                if match:
                    for i in range(6):
                        Stress[i] = match.group(i + 1)
            fid.close()
        return Stress, Energy

    def updateStrain(self, InputStrain, Correct_strain):
        addedstrain = np.zeros([3, 3], "float")
        addedstrain[1, 1] = InputStrain[1]
        addedstrain[2, 2] = InputStrain[2]
#	addedstrain[0,1], addedstrain[1,0] = 0.5 * InputStrain[5] ,0.5 * InputStrain[5]
#	addedstrain[0,2], addedstrain[2,0] = 0.5 * InputStrain[4] ,0.5 * InputStrain[4]
#	addedstrain[1,2], addedstrain[2,1] = 0.5 * InputStrain[3] ,0.5 * InputStrain[3]
        Correct_strain += addedstrain
        Correct_strain = np.mat(Correct_strain)
        return Correct_strain

    def output(self, Stress_original, delta, Correct_strain, Stress):
        with open("DATA", 'a') as ff:
            ff.write("delta %6.4f\t Sxx %10.8f\n" % (delta, Stress_original[0]))
            print("Correct_strain", file=ff)
            print(Correct_strain, file=ff)
            print("Stress original", file=ff)
            print(Stress_original, file=ff)
            print("calculated Stress", file=ff)
            print(Stress, file=ff)
            os.system("cp POSCAR POSCAR%4.3f" % (delta))
            ff.close()
        return

    def gnposcar(self, delta, Correct_strain):
        M = self.getM()
        lattice_constant = self.lattice_constant
        M = self.M
        original_strain = np.matrix(
            [[1 + delta, 0, 0], [0, 1 - 0.1 * delta, 0], [0, 0, 1 + 0.1 * delta]], "float") + Correct_strain
        Base_vector = np.mat('1,0,0;0,1,0;0,0,1', "float")
        Transformed_strain = original_strain
        Transposed_Base = Transformed_strain * Base_vector
        print("Transformed Strain ", original_strain)
        print("transformed Base", Transposed_Base)
        with open("poscar", mode='w') as fout:
            fout.write("bcc\n")
            fout.write("%f\n" % (lattice_constant))
            fout.write(" %22.16f  %22.16f  %22.16f\n" %
                       (Transposed_Base[0, 0], Transposed_Base[0, 1], Transposed_Base[0, 2]))
            fout.write(" %22.16f  %22.16f  %22.16f\n" %
                       (Transposed_Base[1, 0], Transposed_Base[1, 1], Transposed_Base[1, 2]))
            fout.write(" %22.16f  %22.16f  %22.16f\n" %
                       (Transposed_Base[2, 0], Transposed_Base[2, 1], Transposed_Base[2, 2]))
            fout.write("2\n")
            fout.write("Direct\n")
            fout.write("0.0 0.0 0.0\n")
            fout.write("0.5 0.5 0.5\n")
            fout.flush()
            fout.close()
            os.system("mv poscar POSCAR")
        return Transformed_strain


class loop_calculate(object):
    def __init__(self):
        self.execuble = "mpirun vasp"
        self.looptime = 20
        return

    def calculate_tensile(self, Element, lattice_constant, orientation, structure):
        execuble = self.execuble
        #### main function ############
        M = VAtension(Element, lattice_constant, orientation, structure)
        if M.structure == 'bcc':
            for i in range(0, 15):
                Correct_strain = np.mat([[0, 0, 0], [0, 0, 0], [0, 0, 0]], "float")
                delta = i * 0.02
                count = 0
                while True:
                    Transformed_strain = M.gnposcar(delta, Correct_strain)
                    os.system("%s > report.vasp" % (execuble))
                    Stress_original, Energy = M.grabStress()
                    Stress = copy.deepcopy(Stress_original)
                    Stress[0] = 0.0
                    Stress_abs = np.abs(Stress)
                    StressPrime = np.max([Stress_abs[1], Stress_abs[2]])
                    if StressPrime < 1.0:
                        M.output(Stress_original, delta, Correct_strain, Stress)
                        break
                    else:
                        Strain = M.Sij * Stress * 0.06
                        Correct_strain = M.updateStrain(Strain, Correct_strain)
                        with open("moniter.txt", 'a') as fout:
                            print("M.Sij", M.Sij, file=fout)
                            print("loop time",  count, file=fout)
                            print("StressPrime", StressPrime, file=fout)
                            print("original Stress", Stress_original, file=fout)
                            print("Correct_strain", Correct_strain, file=fout)
                            print("transformed Strain", Transformed_strain, file=fout)
                            count += 1
                            fout.close()
        return


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("-s", "--structure",
                      dest="structure",
                      help="fcc bcc or hcp")

    parser.add_option("-e", "--element",
                      dest="Element",
                      help="type of Element",
                      default="Mo")

    parser.add_option("-l",
                      "--lattice_constant",
                      dest="lattice_constant",
                      type="float",
                      help="Lattice constant",
                      default=3.15)

    (options, args) = parser.parse_args()
    Job = loop_calculate()
    Job.calculate_tensile()
