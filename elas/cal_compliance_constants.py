#!/usr/bin/env python
# encoding: utf-8


class cal_para(object):

    def __init__(self):
        return

    def get_cij_sij(self):
        Element_name = []
        C11 = []
        C44 = []
        C12 = []
        Youngs_Modulus = []
        Shear_Modulus = []
        Possion_ratio = []
        A = []
        for i in range(len(Raw) / 9):
            Element_name.append(Raw[9 * i].strip())
            C11.append(float(Raw[(9 * i + 2)]))
            C44.append(float(Raw[(9 * i + 3)]))
            C12.append(float(Raw[(9 * i + 4)]))
            Youngs_Modulus.append(float(Raw[9 * i + 5]))
            Possion_ratio.append(float(Raw[9 * i + 6]))
            Shear_Modulus.append(float(Raw[9 * i + 7]))
            A.append(float(Raw[9 * i + 8]))
        print(Element_name, end=' ')
        print(C11)
        print(C44)
        print(C12)  # Youngs_Modulus, Shear_Modulus, Possion_ratio, A
        c11 = arange(len(Element_name), dtype="float")
        c12 = arange(len(Element_name), dtype="float")
        c44 = arange(len(Element_name), dtype="float")
        s11 = arange(len(Element_name), dtype='float')
        s12 = arange(len(Element_name), dtype='float')
        s44 = arange(len(Element_name), dtype='float')
        for i in range(len(Element_name)):
            # Youngs_Modulus[i] * (1 - Possion_ratio[i]) /
            # (( 1 + Possion_ratio[i]) * (1 - 2 * Possion_ratio[i]))
            c11[i] = C11[i]
            # Youngs_Modulus[i] * Possion_ratio[i] /
            # ((1+Possion_ratio[i])*(1-2*Possion_ratio[i]))
            c12[i] = C12[i]
            c44[i] = C44[i]
            M = zeros([6, 6], "float")
            M[0, 0], M[1, 1], M[2, 2] = c11[i], c11[i], c11[i]
            M[0, 1], M[0, 2], M[1, 0], M[1, 2], M[2, 0], M[
                2, 1] = c12[i], c12[i], c12[i], c12[i], c12[i], c12[i]
            M[3, 3], M[4, 4], M[5, 5] = c44[i], c44[i], c44[i]
            N = inv(M)
            print(N)
            s11[i] = N[0, 0]  # divide(1,float(Youngs_Modulus[i]))
            s12[i] = N[0, 1]  # divide((-Possion_ratio[i]),Youngs_Modulus[i])
            s44[i] = N[4, 4]  # divide(1.,float(Shear_Modulus[i]))
            print(s11[i], s12[i], s44[i])
        fout = open("cij_sij.txt", 'w')
        fout.write(
            "Element  \t c11   \t  c12  \t  c44  \t  s11  \t  s12  \t  s44 \n")
        for i in range(len(Element_name)):
            fout.write("%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n" % (
                Element_name[i], c11[i], c12[i], c44[i], s11[i], s12[i], s44[i]))
        return

    def cal_complience(self, c11, c12, c44):
        import numpy as np
        from numpy.linalg import inv
        M = np.zeros([6, 6], "float")
        M[0, 0], M[1, 1], M[2, 2] = c11, c11, c11
        M[0, 1], M[0, 2], M[1, 0] = c12, c12, c12
        M[1, 2], M[2, 0], M[2, 1] = c12, c12, c12
        M[3, 3], M[4, 4], M[5, 5] = c44, c44, c44
        N = inv(M)
        print(N)
        s11 = N[0, 0]  # divide(1,float(Youngs_Modulus[i]))
        s12 = N[0, 1]  # divide((-Possion_ratio[i]),Youngs_Modulus[i])
        s44 = N[4, 4]  # divide(1.,float(Shear_Modulus[i]))
        print(s11, s12, s44)
        return s11, s12, s44


if __name__ == "__main__":
    from numpy import *
    from numpy.linalg import inv
    from optparse import OptionParser
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-e", "--element", dest="Element",
                      help="type of Element", default="Al")
    parser.add_option("-E", "--Young's modulus",
                      type="float", dest="E", default=70)
    parser.add_option("-v", "--Poissions", dest="Possion",
                      type="float", default=0.35)
    parser.add_option("-s", "--shear", dest="Shear", type="float", default=26)
    parser.add_option("-f", "--file", dest="filename", default="Modulus.txt")
    (options, args) = parser.parse_args()
    filename = options.filename
    fid = open(filename, 'r')
    Raw = fid.readlines()
    fid.close()
