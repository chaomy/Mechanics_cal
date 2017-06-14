#!/usr/bin/env python
# encoding: utf-8

import ase.lattice.cubic as Cubic
import copy, shutil
#from  ase.visualize import view
import ase.io
import numpy as np
import os , re, glob

class cal(object):
    def __init__(self,
                 directions,
                 size,
                 latticeconstant,
                 element,
                 structure,
                 Kpoints):

        self._element = element
        self._directions = directions
        self._lattice_constant = latticeconstant
        self._structure = structure
        self._size = size

        self.exe = "mpirun vasp"

        self._atoms = self.get_SuperCell()
        self._PerfPosition = self._atoms.get_positions()
        self._PerfCells = self._atoms.get_cell()
        self._Kpoints = Kpoints
        return

    def get_SuperCell(self):
        if self._structure == 'bcc':
            atoms = Cubic.BodyCenteredCubic(directions = self._directions,
                                            size = self._size,
                                            latticeconstant = self._lattice_constant,
                                            symbol = self._element,
                                            pbc = (1,1,1))

        elif self._structure == 'fcc':
            atoms = Cubic.FaceCenteredCubic(directions = self._directions,
                                            size = self._size,
                                            latticeconstant = self._lattice_constant,
                                            symbol = self._element,
                                            pbc = (1,1,1))
        """
        ase.io.write(filename="POSCAR",
                     images=Atoms,
                     format='vasp')

        """
        Positions = atoms.get_positions()
        Cells = atoms.get_cell()

        print len(Positions)
        for i in range(8):
            atoms.pop()

        NewPositions = atoms.get_positions()
        NewCells = atoms.get_cell()

        print len(NewPositions)

        print Cells
        print NewCells
        return atoms

    def gndisplacement(self,
                       Cell,
                       Positions,
                       DisplacementVector):

        AtomN = len(Positions)
        Displacement = copy.deepcopy(Positions)

        DisplacementVector = DisplacementVector
        for i in range(AtomN):
            if i <=  9 :
                    Displacement[i] = [0, 0, 0]
            elif i > 9 :
                    Displacement[i] = DisplacementVector
        return Displacement

    def LoopGSF(self):
        RootDir = os.getcwd()
        for j in range(0,1):
            for i in range(0,1):
                DirName = 'Dir-x-%03d-y-%03d'%(i,j)
                os.mkdir(DirName)
                shutil.copy('POTCAR',DirName)

                os.chdir(DirName)
                self.WriteIncar()
                self.WriteKpoints()

                DisplacementVector = [i * 0.0, 0 * j, 0]
                print DisplacementVector
                DisplacementMatrixDirct = self.gndisplacement(self._PerfCells,
                                                              self._PerfPosition,
                                                              DisplacementVector)

                DisplacementMatrix = copy.deepcopy(DisplacementMatrixDirct)

                print "DisplacementMatrix is", DisplacementMatrix

                Localatoms = self._atoms.copy()
                Localatoms.translate(DisplacementMatrix)

                ase.io.write(filename="POSCAR",
                            images=Localatoms,
                            format='vasp')
                #os.system("%s"%(self.exe))
                #(Energy, Stress, Volume) = self.VAgetData()
                #self.output(0.05 * i ,
                        #   Energy,
                        #   Stress)
                os.system("cp POSCAR POSCAR-{0:03d}.vasp".format(i))
                os.chdir(RootDir)
        return

    def genSurface(self):
        self.get_SuperCell()
        self.WriteIncar()
        self.WriteKpoints()
        Localatoms = self._atoms.copy()
        ase.io.write(filename="POSCAR",
                    images=Localatoms,
                    format='vasp')
        os.system("cp POSCAR POSCAR.vasp")
        return

    def runVASP(self):
        DirList = glob.glob("Dir-*")
        print DirList
        RootDir = os.getcwd()
        print len(DirList)
        for i in range(0,1):
            print DirList[i]
            os.chdir(DirList[i])
            os.system("%s"%(self.exe))
            (Energy, Stress, Volume) = self.VAgetData()
            self.output(DirList[i], Energy, Stress)
            os.chdir(RootDir)
        return

    def VAgetData(self):
        Stress=np.arange(6, dtype = "float")
        Stress.shape=([6,1])
        with open("OUTCAR", 'r') as fid:
            real_N = r"(-?\d*\.\d{5})"
            real_Num = r"(-?\d*\.\d*)"
            In = r"\s*"
            for line in fid:
                matchV = \
                re.search(r"\s*volume of cell : \s*" + real_Num,
                         line)
                matchE = \
                re.search(r"\s*energy\s*without\s*entropy=\s*" + \
                         real_Num + In + r"energy\(sigma->0\)\s*=\s*" + \
                         real_Num, line)
                matchS = re.search(r"^\s*in kB\s*" +  real_N + \
                         In + real_N + In + real_N + In + real_N + \
                         In + real_N + In + real_N ,line)
                if matchS :
                    for i in range(6):
                        Stress[i] = float(matchS.group(i + 1))
                if matchE :
                    Energy = float(matchE.group(2))
                if matchV :
                    Volume = float(matchV.group(1))
            fid.close()
        return (Energy, Stress, Volume)

    def WriteKpoints(self):
        with open('KPOINTS', 'w') as fid:
            fid.write("""Automatic mesh
0
Monkhorst-Pack
%d  %d  %d
0   0   0
                    """%(self._Kpoints[0],
                         self._Kpoints[1],
                         self._Kpoints[2]))
        return

    def WriteIncar(self):
        with open('INCAR', 'w') as fid:
            fid.write("""SYSTEM = run.cfg
Start Parameter for This Run:
NWRITE = 2

Electronic Relaxation 1:
PREC   = Accurate
ISYM   = 2 # pay attention !!!
Ionic Relaxation:
NSW    = 2000  ##number of ionic steps taken
NBLOCK = 1
KBLOCK = 1
IBRION = 2
# POTIM  = 0.5
ISIF   = 2 #2 cell shape change no. cell volume change no. 3: allows c
TEBEG  = 300.
SMASS  = 0
LCHARG = .FALSE.
LWAVE = .FALSE.

DOS Related Values:
ISMEAR = 1  #
SIGMA  = 0.4
EMIN   = 10.
EMAX   = 0.
LELF   = T
LVTOT  = T

Electronic Relaxation 2:
#IALGO  = 38
ADDGRID = .TRUE.
LASPH  = .TRUE.
LREAL  = .FALSE.
ISPIN  = 1
ENCUT  = 360.00 # 1.5*(Max ENMAX)
ENAUG  = 560
EDIFF  = 1e-6
LVTOT  = .FALSE.
NPAR   = 4   # close to the square root of # of cores.
""")
            fid.close()
        return

    def output(self,
               Displacement,
               Energy,
               Stress):
        with open("DATA", 'a') as fid:
            fid.write("%s\t%10.8f\t%10.8f\t%10.8f\t"\
                      "%10.8f\t%10.8f\t%10.8f\t%10.8f\n"
                    %(Displacement,Energy,Stress[0],Stress[1],
                      Stress[2],Stress[3],Stress[4],Stress[5]))
            fid.close()
        return

if __name__=='__main__':
    Job = cal(directions =[[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]],
              size = (1,1,18),
              latticeconstant = 3.322404,
              element = 'Nb',
              structure = 'bcc',
              Kpoints = [18,18,1])
    #Job.get_SuperCell()
#    Job.LoopGSF()
#    Job.runVASP()
#    Job.genSurface()
    print Job.VAgetData() 
