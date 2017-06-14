#!/usr/local/bin/python
# encoding: utf-8
#
#
import ase.lattice.compounds as Compounds
import ase.lattice.cubic as Cubic
import copy, shutil, re
import ase.io
import numpy as np
import copy

class Structure(object):
    def __init__(self,
                 directions,
                 size,
                 latticeconstant,
                 element,
                 structure):

        self._element = element
        self._directions = directions
        self._lattice_constant = latticeconstant
        self._structure = structure
        self._size = size
        return

    def gnL12(self):
        #### Pure Cu ###
        la_Cu =  3.49212
        la_Cu3Ag = 3.615000
        atomsCu = Cubic.FaceCenteredCubic(directions = self._directions,
                                        size = (100,60,4),
                                        latticeconstant = \
                                        la_Cu3Ag,
                                        symbol = 'Cu',
                                        pbc = (1,1,1))

        print "Cu is", atomsCu.get_number_of_atoms()
        ase.io.write(filename = "Cu.cfg",
                        images = atomsCu,
                        format = 'cfg')

        atoms = Compounds.AuCu3(directions = self._directions,
                                size = (50,3,8),
                                latticeconstant = \
                                la_Cu3Ag,
                                symbol = self._element,
                                pbc = (1,1,1))

        Cell = atomsCu.get_cell()
        NewCell = copy.deepcopy(Cell)
        #atoms.set_cell(NewCell, scale_atoms = False)
        #atoms.extend(atomsCu)
        # -------------------- add strain --------------------
        Strain = np.mat([[1,   0,   0],
                         [0.0250, 1,0],
                         [0,   0,   1]],'float')
        Positions = atomsCu.get_positions()
        Positions = np.mat(Positions)

        Positions = Positions * Strain

        NewCell = Strain * NewCell

        atomsCu.set_cell(NewCell)
        atomsCu.set_positions(Positions)

        print  "the cell is", atoms.get_cell()

        ase.io.write(filename = "Cu3Au.cfg",
                        images = atomsCu,
                        format = 'cfg')

        ase.io.write(filename = "Cu3Au.xyz",
                        images = atomsCu,
                        format = 'extxyz')
        return
        ##---------- combine --------
#        Cu_Name, Cuxyz, Cell = self.read_cfg('Cu.cfg')
#        Cu3Au_Name, Cu3Auxyz, Cell = self.read_cfg('Cu3Au.cfg')
#        CuNum = len(Cu_Name)
#        Cu3AuNum = len(Cu3Au_Name)
#
#        x1, y1, z1, x2, y2, z2 = [], [], [], [], [], []
#        atom_type1, atom_type2 = [], []
#
#        for i in range(CuNum):
#            x1.append(float(Cuxyz[i][0])*Cell[0]);
#            y1.append(float(Cuxyz[i][1])*Cell[1]);
#            z1.append(float(Cuxyz[i][2])*Cell[2])
#            atom_type1.append(Cu_Name[i])
#
#        for i in range(Cu3AuNum):
#            x2.append(float(Cu3Auxyz[i][0])*Cell[0]);
#            y2.append(float(Cu3Auxyz[i][1])*Cell[1]);
#            z2.append(float(Cu3Auxyz[i][2])*Cell[2])
#            atom_type2.append(Cu3Au_Name[i])
#
#        y2Max = np.max(y2)
#        xnew, ynew, znew, atomtypenew = [], [], [], []
#        for i in range(CuNum):
#            if y1[i] > y2Max:
#                xnew.append(x1[i])
#                ynew.append(y1[i])
#                znew.append(z1[i])
#                atomtypenew.append(Cu_Name[i])
#
#        with open("L1_2.cfg", 'w') as fid:
#            fid.write("Number of particles = %d \n" %(len(xnew) +
#                Cu3AuNum))
#            fid.write("A = 1.0 Angstrom\n")
#            fid.write("""H0(1,1) = %f A
#H0(1,2) = 0.000000 A
#H0(1,3) = 0.000000 A
#H0(2,1) = 0.000000 A
#H0(2,2) = %f A
#H0(2,3) = 0.000000 A
#H0(3,1) = 0.000000 A
#H0(3,2) = 0.000000 A
#H0(3,3) = %f A
#.NO_VELOCITY.
#entry_count = 3
#"""%(Cell[0], Cell[1], Cell[2]))
#            for i in range(Cu3AuNum):
#                if atom_type2[i] == 'Cu':
#                    fid.write("63.546000\n")
#                    fid.write("Cu\n")
#                    fid.write("%f  %f  %f\n"
#                             %(x2[i]/Cell[0], y2[i]/Cell[1],
#                                 z2[i]/Cell[2]))
#                elif atom_type2[i] == 'Au':
#                    fid.write("196.9665\n")
#                    fid.write("Ag\n")
#                    fid.write("%f  %f  %f\n"
#                             %(x2[i]/Cell[0], y2[i]/Cell[1],
#                                 z2[i]/Cell[2]))
#            for i in range(len(xnew)):
#                    fid.write("63.546000\n")
#                    fid.write("Cu\n")
#                    fid.write("%f  %f  %f\n"
#                             %(xnew[i]/Cell[0], ynew[i]/Cell[1],
#                                 znew[i]/Cell[2]))

    def read_cfg(self, filename):
        with open(filename, 'r') as fid:
            Raw = fid.read()
        print Raw
        xyz = re.compile(r"([+\-]?\d*\.\d*e[+\-]?\d+) ([+\-]?\d*\.\d*e[+\-]?\d+) ([+\-]?\d*\.\d*e[+\-]?\d+)")
        element =  re.compile(r"(Cu|Au)")
        xhi = re.compile(r"H0\(1,1\)\s*=\s*(\d*.\d*)")
        yhi = re.compile(r"H0\(2,2\)\s*=\s*(\d*.\d*)")
        zhi = re.compile(r"H0\(3,3\)\s*=\s*(\d*.\d*)")

        xyzdata = xyz.findall(Raw)
        elementdata = element.findall(Raw)
        Cell = []
        Cell.append(float(xhi.findall(Raw)[0]))
        Cell.append(float(yhi.findall(Raw)[0]))
        Cell.append(float(zhi.findall(Raw)[0]))
        print Cell
        print float(xyzdata[0][0])
        print float(xyzdata[1][0])
        return (elementdata, xyzdata, Cell)


if __name__=='__main__':
    Job = Structure(directions =[[-1, 1, 0],
                                 [1, 1, 0],
                                 [0, 0, -1]],
              size = (2,2,2),
              latticeconstant = 3.690900,
              element = ('Au', 'Cu'),
              structure = 'L12')
    Job.gnL12()

