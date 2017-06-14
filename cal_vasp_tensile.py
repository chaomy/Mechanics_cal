#!/usr/local/bin/python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_vasp_tensile.py
#
###################################################################
#
# Purpose :  multithreads to calculate vasp tensile
#
# Creation Date :
# Last Modified : Thu Mar 30 23:56:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import glob
import os
import numpy as np
import shutil
import sys
import ase
import matplotlib.pylab as plt

try:
    import gn_config
    import get_data
    import gn_incar
    import gn_kpoints
    import gn_pbs

except ImportError:
    print("error during import")


class vasp_itensile(gn_config.bcc,
                    get_data.get_data,
                    gn_incar.gn_incar,
                    gn_pbs.gn_pbs,
                    gn_kpoints.gn_kpoints):

    def __init__(self):
        gn_incar.gn_incar.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gn_config.bcc.__init__(self)

        self.root_dir = os.getcwd()
        self.set_element('Nb')
        self.num = 21
        return

    def prepare_dislocation_vasp_infiles(self, dirname=None):
        self.set_incar_type("dft")

        self.write_incar()
        #  self.set_diff_kpoints([18, 18, 18]);  # tpath
        #  self.set_diff_kpoints([18, 14, 14]);  # opath

        self.set_diff_kpoints([9, 9, 9])

        self.write_kpoints()

        self.set_nnodes(1)
        self.set_ppn(8)
        self.set_pbs_type('va')
        self.set_wall_time(90)

        if dirname is None:
            self.set_job_title("tpath")
        else:
            self.set_job_title(dirname)

        self.set_main_job("mpirun vasp")
        self.write_pbs()

        os.system("cp ../POTCAR  .")
        return

    def recal_relaxed_structure(self):
        for i in range(0, 26):

            file = "POSCAR%.3f" % (0.01 * i)
            dirname = "dir-%03d" % (i)

            if (not os.path.isdir(dirname)):
                os.mkdir(dirname)

            shutil.copy(file, dirname)
            #  shutil.copy("POTCAR", dirname);

            os.chdir(dirname)

            shutil.copy(file, "POSCAR")
            shutil.copy(file,  self.root_dir)

            self.prepare_dislocation_vasp_infiles(dirname)

            os.chdir(self.root_dir)
        return

    def get_energy(self):
        for i in range(26):
            dirname = "dir-%03d" % (i)
            os.chdir(dirname)
            #  self.vasp_energy_stress_vol_quick();

            (energy, stress, vol) = self.vasp_energy_stress_vol()

            outstr = "%.6f %.6f %.6f %.6f %.6f %.6f %.6f" % (energy,
                                                             stress[0], stress[1], stress[2],
                                                             stress[3], stress[4], stress[5])
            print outstr

            os.chdir(self.root_dir)

        return

    def loopFitStrain(self):
        dir_list = glob.glob("dir-*")

        for i in range(len(dir_list)):
            work_dir = dir_list[i]
            os.chdir(work_dir)

            if sys.argv[1] == "sub":
                os.system("qsub va.pbs")
            else:
                self.prepare_dislocation_vasp_infiles(work_dir)

            os.chdir(self.root_dir)
        return

    def read_Strain_stress(self):
        (energy, stress, vol) = self.vasp_energy_stress_vol()
        #  (atom_number, base, comment, atom_position) = \
        #  self.read_vasp_poscar("CONTCAR");
        atoms = ase.io.read("CONTCAR", format == "vasp")

        #  print ("energy", energy);
        #  print ("stress", stress);
        return (energy, stress, atoms.get_cell())

    def relative_to_cell(self, rel_pos, base):
        pos = [0, 0, 0]

        pos[0] = rel_pos[0] * base[0, 0] + rel_pos[1] * base[0, 1] + rel_pos[2] * base[0, 2]
        pos[1] = rel_pos[0] * base[1, 0] + rel_pos[1] * base[1, 1] + rel_pos[2] * base[1, 2]
        pos[2] = rel_pos[0] * base[2, 0] + rel_pos[1] * base[2, 1] + rel_pos[2] * base[2, 2]

        return pos

    def pre_dft_iten_data(self, strain=0.02, tag="tpath"):
        ##### we expand the box to be 2 by 2 by 2 ####

        msize = [2, 2, 2]
        #  atom_num = 4 * msize[0] * msize[1] * msize[2];

        (a, b, c) = [1. / msize[i] for i in range(3)]
        prelist = [2, 6, 10, 14, 18]

        #  if tag == "opath":
        for j in prelist:
            dirname = "dir-%.3f-%.3d" % (strain, j)
            dirnamecnt = "cnt-%.3f-%.3d" % (strain, j)

            if (not os.path.isdir(dirnamecnt)):
                os.mkdir(dirnamecnt)

            ### the returned base already has strain ###
            os.chdir(dirname)
            (energy, stress, base) = self.read_Strain_stress()

            ### exapnd the box ###
            for i in range(3):
                base[i] = base[i] * msize[i]

            ### assign atoms positions ###
            temp_pos = [0, 0, 0]

            ### initialize atoms ###
            atoms = ase.Atoms(cell=base, pbc=[1, 1, 1])

            if tag == "opath":
                for z in range(msize[0]):
                    for y in range(msize[1]):
                        for x in range(msize[2]):
                            temp_pos[0] = x * a
                            temp_pos[1] = y * b
                            temp_pos[2] = z * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

                            temp_pos[0] = (x + 0.5) * a
                            temp_pos[1] = (y + 0.5) * b
                            temp_pos[2] = z * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

                            temp_pos[0] = x * a
                            temp_pos[1] = (y + 0.5) * b
                            temp_pos[2] = (z + 0.5) * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

                            temp_pos[0] = (x + 0.5) * a
                            temp_pos[1] = y * b
                            temp_pos[2] = (z + 0.5) * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

            if tag == "tpath":
                for z in range(msize[0]):
                    for y in range(msize[1]):
                        for x in range(msize[2]):

                            temp_pos[0] = x * a
                            temp_pos[1] = y * b
                            temp_pos[2] = z * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

                            temp_pos[0] = (x + 0.5) * a
                            temp_pos[1] = (y + 0.5) * b
                            temp_pos[2] = (z + 0.5) * c
                            pos = self.relative_to_cell(temp_pos, base)
                            atoms.append(ase.Atom("Nb", pos))

                os.chdir(self.root_dir)
                os.chdir(dirnamecnt)

                ase.io.write("POSCAR", images=atoms, format='vasp')
                os.system("cp POSCAR POSCAR.vasp")
                self.prepare_dislocation_vasp_infiles()

                os.chdir(self.root_dir)
        return

    def loopCollectStrain(self, strain=None):
        num = self.num

        #  for i in range():
        strn_listT = np.zeros(num)
        strn_listO = np.zeros(num)

        engy_listT = np.zeros(num)
        engy_listO = np.zeros(num)

        sxx_listT = np.zeros(num)
        sxx_listO = np.zeros(num)

        strss_listT = []
        strss_listO = []
        base_listT = []
        base_listO = []

        if strain is None:
            strain = 0.04

        lattice = 3.322404
        sqrt2_lattice = np.sqrt(2) * lattice

        for j in range(21):
            ###############  TP  ##############
            dirnameT = "FitStrainsTP/dir-%.3f-%.3d" % (strain, j)

            os.chdir(dirnameT)
            (energy, stress, base) = self.read_Strain_stress()
            print base

            strn_listT[j] = base[1, 1]
            engy_listT[j] = energy
            sxx_listT[j] = stress[0]

            strss_listT.append(stress)
            base = base / lattice
            base_listT.append(base)
            os.chdir(self.root_dir)

            ###############  OP  ##############
            dirnameO = "FitStrainsOP/dir-%.3f-%.3d" % (strain, j)
            os.chdir(dirnameO)
            (energy, stress, base) = self.read_Strain_stress()

            strn_listO[j] = base[1, 1]
            engy_listO[j] = energy
            sxx_listO[j] = stress[0]

            strss_listO.append(stress)
            base[0, :] = base[0, :] / lattice;
            base[1, :] = base[1, :] / sqrt2_lattice;
            base[2, :] = base[2, :] / sqrt2_lattice;

            base_listO.append(base)
            os.chdir(self.root_dir)

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(strn_listO, engy_listO)

        ax2 = fig.add_subplot(212)
        ax2.plot(strn_listO, sxx_listO)

        plt.savefig("strain-%.3f.png" % (strain))

        ####### be careful the units in VASP is kbar ######
        with open("log-%.3f" % (strain), "w") as fid:
            fid.write("    if (strain == %.3f){\n" % (strain))

            for i in range(num):

                fid.write("""
        strcpy(temp.name, "Nb");
        temp.tpathStrainM[0][0]= %10.8f; temp.tpathStrainM[0][1]= %10.8f; temp.tpathStrainM[0][2]= %10.8f;
        temp.tpathStrainM[1][0]= %10.8f; temp.tpathStrainM[1][1]= %10.8f; temp.tpathStrainM[1][2]= %10.8f;
        temp.tpathStrainM[2][0]= %10.8f; temp.tpathStrainM[2][1]= %10.8f; temp.tpathStrainM[2][2]= %10.8f;

        temp.opathStrainM[0][0]= %10.8f; temp.opathStrainM[0][1]= %10.8f; temp.opathStrainM[0][2]= %10.8f;
        temp.opathStrainM[1][0]= %10.8f; temp.opathStrainM[1][1]= %10.8f; temp.opathStrainM[1][2]= %10.8f;
        temp.opathStrainM[2][0]= %10.8f; temp.opathStrainM[2][1]= %10.8f; temp.opathStrainM[2][2]= %10.8f;

        temp.tpathEnergy = %10.8f;
        temp.tpathStressV[0]=  %10.8f; temp.tpathStressV[1] =  %10.8f; temp.tpathStressV[2]=  %10.8f;
        temp.tpathStressV[3]=  %10.8f; temp.tpathStressV[4] =  %10.8f; temp.tpathStressV[5]=  %10.8f;

        temp.opathEnergy = %10.8f;
        temp.opathStressV[0]=  %10.8f; temp.opathStressV[1] =  %10.8f; temp.opathStressV[2]=  %10.8f;
        temp.opathStressV[3]=  %10.8f; temp.opathStressV[4] =  %10.8f; temp.opathStressV[5]=  %10.8f;

        genItenlist.push_back(temp);
        """ % (
                    base_listT[i][0, 0], base_listT[i][0, 1], base_listT[i][0, 2],
                    base_listT[i][1, 0], base_listT[i][1, 1], base_listT[i][1, 2],
                    base_listT[i][2, 0], base_listT[i][2, 1], base_listT[i][2, 2],

                    base_listO[i][0, 0], base_listO[i][0, 1], base_listO[i][0, 2],
                    base_listO[i][1, 0], base_listO[i][1, 1], base_listO[i][1, 2],
                    base_listO[i][2, 0], base_listO[i][2, 1], base_listO[i][2, 2],

                    engy_listT[i],
                    strss_listT[i][0], strss_listT[i][1], strss_listT[i][2],
                    strss_listT[i][3], strss_listT[i][4], strss_listT[i][5],

                    engy_listO[i],
                    strss_listO[i][0], strss_listO[i][1], strss_listO[i][2],
                    strss_listO[i][3], strss_listO[i][4], strss_listO[i][5]))
            fid.write("\n")
            fid.write("    }\n")
            fid.close()
        return

    def looploopCollectStrain(self):
        tag = "tpath"
        for i in range(12):
            strain = i * 0.02
            #  self.loopCollectStrain(strain);
            self.pre_dft_iten_data(strain, tag)
            #  filename = "log-%.3f"%(strain);
            #  os.system("cat %s >> lmpStrainDataDFT.cpp"%(filename));
        return

    def sub_jobs(self):
        dir_list = glob.glob("cnt-*")
        for mdir in dir_list:
            os.chdir(mdir)

            os.system("qsub va.pbs")

            os.chdir(self.root_dir)

        return


if __name__ == '__main__':
    job = vasp_itensile()
    #  job.recal_relaxed_structure();
    #  job.get_energy();

    if ((sys.argv[1] == "pre") or (sys.argv[1] == "sub")):
        job.loopFitStrain()

    if sys.argv[1] == "pos":
        job.loopCollectStrain()

    if sys.argv[1] == "loop":
        job.looploopCollectStrain()

    if sys.argv[1] == "pre_iten":
        job.pre_dft_iten_data()

    if sys.argv[1] == "sub":
        job.sub_jobs()
