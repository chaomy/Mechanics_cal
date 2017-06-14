#!/usr/local/bin/python
# encoding: utf-8
#
import copy
import os
import numpy as np
from multiprocessing import Pool
import shutil

try:
    import  get_data
    import  output_data
    import  gn_lmp_infile

except ImportError:
    print "error during import"

__version__ = 0.01
__author__ = 'Chaoming Yang'
__all__=['md_tensile']

def unwrap_self_run_lammps(arg, **kwarg):
    return  md_loop_tensile.lammps_job(*arg, **kwarg)

class md_tensile(get_data.get_data,
                 output_data.output_data,
                 gn_lmp_infile.gn_md_infile):

    def __init__(self,element,
                 lattice_constant,
                 size,
                 orientation,structure):

        get_data.get_data.__init__(self)
        output_data.output_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)

        self.element = element
        self.lattice_constant = lattice_constant
        self.size = size
        self.orientation = orientation
        self.M = self.get_M()
        self.Sij = self.get_Sij()
        self.structure = structure

        self.lx,self.ly,self.lz = 0, 0, 0
        self.addedstrain = np.mat([[0,0,0],
                                   [0,0,0],
                                   [0,0,0]],"float")
        self.stress_original = np.array([0,0,0,0,0,0],
                                          dtype="float")
        return

    def update_strain(self,
                     InputStrain,
                     Correct_strain):
        addedstrain = copy.deepcopy(self.addedstrain)
        addedstrain[1,1] = InputStrain[1]
        addedstrain[2,2] = InputStrain[2]
        Correct_strain += addedstrain

        return Correct_strain

    def gn_bcc_tpath(self, delta, Correct_strain):
        element = self.element
        M = self.get_M()
        lattice_constant = self.lattice_constant
        size = self.size
        M = self.M
        original_strain=np.matrix([[1+delta,0,0],
                                   [0,1+0.2*delta,0],
                                   [0,0,1-0.2*delta]],
                                   "float") + Correct_strain
        Transformed_strain= M.transpose() * original_strain * M
        Base_vector = np.matrix([[1,0,0],
                                 [0,1,0],
                                 [0,0,1]])
        Transposed_Base =  Transformed_strain * Base_vector

        xsize=size[0]
        ysize=size[1]
        zsize=size[2]
        a = 1./ xsize
        b = 1./ ysize
        c = 1./ zsize
        Lengthx = Transposed_Base[0,:] * Transposed_Base[0,:].transpose()
        Lengthy = Transposed_Base[1,:] * Transposed_Base[1,:].transpose()
        Lengthz = Transposed_Base[2,:] * Transposed_Base[2,:].transpose()
        Lengthx = np.sqrt(Lengthx) ; Lengthx *= lattice_constant
        Lengthy = np.sqrt(Lengthy) ; Lengthy *= lattice_constant
        Lengthz = np.sqrt(Lengthz) ; Lengthz *= lattice_constant

        xlo, xhi = 0.0, Lengthx * xsize
        ylo, yhi = 0.0, Lengthy * ysize
        zlo, zhi = 0.0, Lengthz * zsize

        AtomNumber = int(xsize * ysize * zsize) * 2
        XDirection, YDirection , ZDirection = [], [] ,[]

        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XDirection.append((x) * a * Lengthx)
                    YDirection.append((y) * b * Lengthy)
                    ZDirection.append((z) * c * Lengthz)

                    XDirection.append((x + 0.5) * a * Lengthx)
                    YDirection.append((y + 0.5) * b * Lengthy)
                    ZDirection.append((z + 0.5) * c * Lengthz)

        filename = "lattice.txt" #%(element,lattice_constant)
        with open(filename,mode="w") as fout:
            fout.write("#bcc_structure_input for lammps\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0, 0, 0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                        %(i+1,
                         XDirection[i],
                         YDirection[i],
                         ZDirection[i]))
        fout.close()
        os.system("cp lattice.txt lattice_%5.3f.txt"%(delta))
        ########### generate cfg #########################
        XXDirection , YYDirection , ZZDirection = [],[],[]
        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XXDirection.append((x) * a )
                    YYDirection.append((y) * b )
                    ZZDirection.append ((z) * c)

                    XXDirection.append((x + 0.5) * a )
                    YYDirection.append((y + 0.5) * b )
                    ZZDirection.append((z + 0.5) * c )

        curdir=os.getcwd()
        os.chdir("mycfg")
        filename = "cfg_%5.4f.cfg"%(delta)
        with open(filename,mode="w") as fout:
            fout.write("Number of particles = %d\n"%(AtomNumber))
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
    """%(Transposed_Base[0,0]*xsize*Lengthx,
         Transposed_Base[1,1]*ysize*Lengthy,
         Transposed_Base[2,2]*zsize*Lengthz,
         Transposed_Base[0,1]*xsize,
         Transposed_Base[0,2]*xsize,
         Transposed_Base[1,0]*xsize,
         Transposed_Base[1,2]*xsize,
         Transposed_Base[2,0]*xsize,
         Transposed_Base[2,1]*xsize))

            for i in range(AtomNumber):
                fout.write("95.94\t%s\t%7.7f\t %7.7f\t %7.7f\t0\t0\t0\n"
                        %(element,
                         XXDirection[i],
                         YYDirection[i],
                         ZZDirection[i]))
        fout.close()
        os.chdir(curdir)
        return

    def gn_bcc_opath(self,delta,Correct_strain):
        element = self.element
        M = self.get_M()
        lattice_constant = self.lattice_constant
        size = self.size
        M = self.M
        original_strain=np.matrix([[1 + delta, 0, 0],
                                  [0,1 + 1.2*delta,0],
                                  [0,0,1 - 1.2*delta]],
                                  "float") + Correct_strain
        Transformed_strain= M.transpose() * original_strain * M
        Base_vector = np.matrix([[1,        0,          0],
                                [ 0,np.sqrt(2), 	0],
                                [ 0,        0,np.sqrt(2)]])
        Transposed_Base =  Transformed_strain * Base_vector
        xsize=size[0]
        ysize=size[1]
        zsize=size[2]
        a = 1./ xsize
        b = 1./ ysize
        c = 1./ zsize
        Lengthx = Transposed_Base[0,:] * Transposed_Base[0,:].transpose()
        Lengthy = Transposed_Base[1,:] * Transposed_Base[1,:].transpose()
        Lengthz = Transposed_Base[2,:] * Transposed_Base[2,:].transpose()
        Lengthx = np.sqrt(Lengthx) ; Lengthx *= lattice_constant
        Lengthy = np.sqrt(Lengthy) ; Lengthy *= lattice_constant
        Lengthz = np.sqrt(Lengthz) ; Lengthz *= lattice_constant

        xlo, xhi = 0.0, Lengthx * xsize
        ylo, yhi = 0.0, Lengthy * ysize
        zlo, zhi = 0.0, Lengthz * zsize

        AtomNumber = int(xsize * ysize * zsize) * 4
        XDirection, YDirection , ZDirection = [], [] ,[]

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
                    ZDirection.append((z + 0.5) * c * Lengthz )

                    XDirection.append((x + 0.5) * a * Lengthx)
                    YDirection.append((y + 0.5) * b * Lengthy)
                    ZDirection.append((z) * c * Lengthz)

        filename = "lattice.txt" #%(element,lattice_constant)
        with open(filename,mode="w") as fout:
            fout.write("#bcc_structure_input for lammps\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                        %(i+1,
                         XDirection[i],
                         YDirection[i],
                         ZDirection[i]))
        fout.close()
        os.system("cp lattice.txt lattice_%5.3f.txt"%(delta))
        ########### generate cfg #########################
        XXDirection , YYDirection , ZZDirection = [],[],[]
        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    XXDirection.append((x) * a )
                    YYDirection.append((y) * b )
                    ZZDirection.append((z) * c)

                    XXDirection.append((x) * a)
                    YYDirection.append((y + 0.5) * b )
                    ZZDirection.append((z + 0.5) * c )

                    XXDirection.append((x + 0.5) * a )
                    YYDirection.append((y) * b )
                    ZZDirection.append((z + 0.5) * c )

                    XXDirection.append((x + 0.5) * a )
                    YYDirection.append((y + 0.5) * b )
                    ZZDirection.append ((z ) * c )
        curdir=os.getcwd()
        os.chdir("mycfg")
        filename = "cfg_%5.4f.cfg"%(delta)
        with open(filename,mode="w") as fout:
            fout.write("Number of particles = %d\n"%(AtomNumber))
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
"""%(Transposed_Base[0,0]*xsize*Lengthx,
     Transposed_Base[1,1]*ysize*Lengthy,
     Transposed_Base[2,2]*zsize*Lengthz,
     Transposed_Base[0,1]*xsize,
     Transposed_Base[0,2]*xsize,
     Transposed_Base[1,0]*xsize,
     Transposed_Base[1,2]*xsize,
     Transposed_Base[2,0]*xsize,
     Transposed_Base[2,1]*xsize))
            for i in range(AtomNumber):
                fout.write("95.94\t%s\t%7.7f\t %7.7f\t %7.7f\t0\t0\t0\n"
                        %(element,
                        XXDirection[i],
                        YYDirection[i],
                        ZZDirection[i]))
        fout.close()
        os.chdir(curdir)
        return

class md_loop_tensile(md_tensile):
    def __init__(self,
                 element = None,
                 lattice_constant = None,
                 size = None,
                 orientation = None,
                 structure = None,
                 in_potential = None):

        self._flux_exe =  'lmp_mpi < in.stat_tensile'

        self._looptime = 11
        self._increment = 0.02
        self._Stress_ThrValue = 1.0

        if in_potential is not None:
            self._tensile_potential = in_potential
        else:
            self._tensile_potential = 'dummy.lammps.EAM'  # default

        if element is not None:
            self._element = element
        else:
            self._element = 'Nb'

        if lattice_constant is not None:
            self._lattice_constant = lattice_constant
        else:
            self._lattice_constant = 3.30

        self._size = size

        if orientation is not None:
            self._orientation = orientation
        else:
            self._orientation = [[1, 0, 0],
                                 [0, 1, 0],
                                 [0, 0, 1]]

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

    def get_Strain(self,
                   Sij,
                   Stress,
                   StressPrime):
        StressPrime = float(StressPrime)
        Coeff = 0.06 + 0.06 * StressPrime ** 0.20 + \
                         0.03 * StressPrime ** 0.19 + \
                         0.01 * StressPrime ** 0.1
        if StressPrime > 20 :
            Coeff =   0.08 + 0.01 * StressPrime ** 0.19 + \
                             0.09 * StressPrime ** 0.1
        Strain = -Sij * Stress * Coeff
        return Strain

    def calculate_tensile(self,
                          option):
        if option == 'TP' or option == 'tpath':
            Directory="TP_%s_%s_%d"%(self._element,
                                     self._lattice_constant,
                                     self._size[1])
            output_data_file = os.path.join(self.root_dir, 'TP-DATA')

        elif option == 'OP' or option == 'opath':
            Directory="OP_%s_%s_%d"%(self._element,
                                     self._lattice_constant,
                                     self._size[1])
            output_data_file = os.path.join(self.root_dir, 'OP-DATA')

        if os.path.isdir(Directory):
            shutil.rmtree(Directory)

        os.mkdir(Directory)
        shutil.copy("./in.stat_tensile",Directory)
        shutil.copy(self._tensile_potential, Directory)
        os.system("cp ~/My_cal/Generate_config/cij_sij.txt  %s"%(Directory))
        os.chdir(Directory)
        os.mkdir("cfg")
        os.mkdir("mycfg")
        os.mkdir("dump")
        os.mkdir("restart")
        os.mkdir("xyz")

        #### main function ############
        Correct_strain = np.mat([[0,0,0],
                                 [0,0,0],
                                 [0,0,0]],"float")

        md_tensile.__init__(self,
                            self._element,
                            self._lattice_constant,
                            self._size,
                            self._orientation,
                            self._structure)

        for i in range(self._looptime):
            delta = i * self._increment
            count = 0
            while True :
                if option == 'TP':
                    self.gn_bcc_tpath(delta, Correct_strain)
                elif option == 'OP':
                    self.gn_bcc_opath(delta, Correct_strain)

                os.system("cat lattice.txt >> backup.txt")
                os.system("%s > Log_MD"%(self._flux_exe))

                self.lx, self.ly, self.lz = self.md_get_lx_ly_lz()
                self.stress_original = self.md_get_stress()
                Stress = copy.deepcopy(self.stress_original)

                Stress[0] = 0.0
                Stress_abs = np.abs(Stress)
                StressPrime = np.max(Stress_abs)

                if StressPrime < self._Stress_ThrValue:
                    self.output_md_tensile(delta)
                    break
                elif count > 500:
                    self.output_md_tensile(delta)
                    break
                else :
                    Strain = self.get_Strain(self.Sij,
                                             Stress,
                                             StressPrime)
                    Correct_strain = self.update_strain(Strain,
                                                    Correct_strain)
                    with open('monitor.txt', 'a') as fid:
                        print >> fid, "delta ", delta
                        print >> fid, "Run times", count
                        print >> fid, "stress_original", self.stress_original
                        print >> fid, "StressPrime ", StressPrime
                        print >> fid, "Correct_strain", Correct_strain
                        fid.close()
                    count += 1

        shutil.copy('DATA',output_data_file)
        os.chdir(self.root_dir)
        return

    def lammps_job(self,
                  option,
                  job):

        if option == 'TP' or option == 'tpath':
            Directory="TP_%s_%d"%(self._element,
                                  job)
            output_data_file = os.path.join(self.root_dir, 'TP-DATA')
        elif option == 'OP' or option == 'opath':
            Directory="OP_%s_%d"%(self._element,
                                  job)
            output_data_file = os.path.join(self.root_dir, 'OP-DATA')

        if os.path.isdir(Directory):
            shutil.rmtree(Directory)

        os.mkdir(Directory)
        #shutil.copy("./in.stat_tensile", Directory)
        shutil.copy(self._tensile_potential, Directory)
        os.system("cp ~/My_cal/Generate_config/cij_sij.txt %s"%(Directory))
        os.chdir(Directory)

        self.gn_md_tensile(potential_file = self._tensile_potential,
                           element = self._element)
        os.mkdir("cfg")
        os.mkdir("mycfg")
        os.mkdir("dump")
        os.mkdir("restart")
        os.mkdir("xyz")

        if os.path.isfile(output_data_file):
            os.system(": > %s"%(output_data_file))
        #### main function ############
        Correct_strain = np.mat([[0,0,0],
                                [0,0,0],
                                [0,0,0]],"float")

        md_tensile.__init__(self,
                            self._element,
                            self._lattice_constant,
                            self._size,
                            self._orientation,
                            self._structure)

        delta = job * 0.02
        count = 0
        while True :
            if option == 'TP':
                self.gn_bcc_tpath(delta, Correct_strain)
            elif option == 'OP':
                self.gn_bcc_opath(delta, Correct_strain)

            os.system("cat lattice.txt >> backup.txt")
            os.system("%s > Log_MD"%(self._flux_exe))

            self.lx, self.ly, self.lz = self.md_get_lx_ly_lz()
            self.stress_original = self.md_get_stress()
            Stress = copy.deepcopy(self.stress_original)

            Stress[0] = 0.0
            Stress_abs = np.abs(Stress)
            StressPrime = np.max(Stress_abs)

            if StressPrime < self._Stress_ThrValue:
                self.output_md_tensile(delta)
                break
            elif count > 500:
                self.output_md_tensile(delta)
                break
            else :
                Strain = self.get_Strain(self.Sij,
                                         Stress,
                                         StressPrime)
                Correct_strain = self.update_strain(Strain,
                                                    Correct_strain)
                with open('monitor.txt', 'a') as fid:
                    print >> fid, "delta ", delta
                    print >> fid, "Run times", count
                    print >> fid, "stress_original", self.stress_original
                    print >> fid, "StressPrime ", StressPrime
                    print >> fid, "Correct_strain", Correct_strain
                    fid.close()
                count += 1
        os.system("cat DATA >> %s"%(output_data_file))
        os.chdir(self.root_dir)
        return (delta, self.stress_original[0,0])

    def Mpi_tensile(self,
                    options):
        pool = Pool(processes = self._looptime)
        List = np.arange(self._looptime)
        Results = pool.map(unwrap_self_run_lammps,
                          zip([self] * len(List),
                              [options] * len(List),
                              List))
        Strain = []
        Stress = []
        for i in range(self._looptime):
            Strain.append(Results[i][0])
            Stress.append(Results[i][1])
        return (Strain, Stress)

if __name__ == "__main__":
    M = md_loop_tensile(element = 'Nb',
                        lattice_constant = 3.307,
                        size = [1,1,1],
                        orientation = [1,0,0,0,1,0],
                        structure = 'bcc',
                        in_potential = "Nb.eam.alloy.webarchive")

    (Strain , Stress) =  M.Mpi_tensile('TP')
    (Strain2 , Stress2) = M.Mpi_tensile('OP')
