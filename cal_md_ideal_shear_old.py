#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_md_ideal_shear_old.py
# 
###################################################################
#
# Purpose :
# 
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

try:
    import

except ImportError:
    print("error during import")


    def md_ideal_shear(self, tag='clc', opt='lmp'):
        basis = self.basis
        npts = 15
        data = np.ndarray([npts, 8])
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)

            delta = 0.02 * i
            strain = np.mat([[1.0,   0.0,  0.0],
                             [-delta, 1.0,  0.0],
                             [0.0,   0.0,  1.0]])

            #  tensile = np.mat([[1.0 + delta,   0.0, 0.0],
            #  [0.0,   1.0 - 0.1 * delta, 0.0],
            #  [0.0,   0.0,  1.0 - 0.1 * delta]])

            if tag == 'run':
                if opt == 'vasp':
                    self.alat = self.lat_vasp
                elif opt == 'lmp':
                    self.alat = self.lat_lmp

                new_strain = basis.transpose() * strain * basis

                self.gn_primitive_lmps(new_strain)
                self.gn_convention(new_strain)

                if opt == 'vasp':
                    # vasp
                    os.system("cp POSCAR_a  {}/POSCAR".format(dirname))
                    os.system("mv POSCAR_a  POSCAR_{:03d}".format(i))
                    os.system("cp KPOINTS   {}".format(dirname))
                    os.system("cp INCAR   {}".format(dirname))
                    os.system("cp POTCAR  {}".format(dirname))
                    os.system("cp va.pbs  {}".format(dirname))

                elif opt == 'lmp':
                    # lmp
                    os.system("lmp_mpi -i  in.init")
                    os.system("mv init.txt  init_{:03d}".format(i))
                    os.system("mv out.txt {}".format(dirname))

            elif tag == 'clc':
                os.chdir(dirname)
                if opt == 'vasp':
                    (engy, stress, vol) = self.vasp_energy_stress_vol()
                    stress *= 0.1
                    stress = self.trans_stress_to_cartesian(stress)
                    data[i][0] = delta
                    data[i][1] = engy
                    data[i][2:] = stress[:, 0]

                elif opt == 'lmp':
                    raw = np.loadtxt("out.txt")
                    engy = raw[0]
                    stress = raw[1:]
                    #  stress = self.trans_stress_to_cartesian(stress)
                    data[i][0] = delta
                    data[i][1] = engy
                    data[i][2:] = stress[:]

                os.chdir(self.root)
        data[:, 1] = data[:, 1] - data[0, 1]
        print data
        return

    def trans_stress_tensile(self):
        delta = 0.02 * 3
        lmpmtx = np.mat(np.zeros([3, 3]))

        #  tensile = np.mat([[1.0 + delta,   0.0, 0.0],
        #  [0.0,   0.97346549918, 0.0],
        #  [0.0,   0.0, 0.97667253459]])

        stsslmp = np.array([-0.48763, -6.18858, -6.023485, 0.94804, -2.579288, 1.057818])

        # Tij = aik * Tkl * ajl
        lmpmtx[0, 0] = stsslmp[0]
        lmpmtx[1, 1] = stsslmp[1]
        lmpmtx[2, 2] = stsslmp[2]

        lmpmtx[1, 0] = stsslmp[3]
        lmpmtx[2, 1] = stsslmp[4]
        lmpmtx[2, 0] = stsslmp[5]

        lmpmtx[0, 1] = stsslmp[3]
        lmpmtx[1, 2] = stsslmp[4]
        lmpmtx[0, 2] = stsslmp[5]

        invbas = np.linalg.inv(bas)
        print invbas * lmpmtx * invbas.transpose()
        return

    def trans_stress(self):
        basis = self.basis

        lmpmtx = np.mat(np.zeros([3, 3]))
        vaspmtx = np.mat(np.zeros([3, 3]))
        conmtx = np.mat(np.zeros([3, 3]))

        ###################################################################
        # lmps
        ###################################################################
        lmpcell = np.mat([[2.9006,  0.,      0.],
                          [-0.905,   2.7558,  0.],
                          [-0.977,  -1.3491,  2.294]])
        invlmp = np.linalg.inv(lmpcell)
        invlmp /= np.linalg.norm(invlmp)
        stsslmp = np.array([-1.77674, -1.77674, 1.070929, -0.389782, 0.059423, 0.059423])

        lmpmtx[0, 0] = stsslmp[0]
        lmpmtx[1, 1] = stsslmp[1]
        lmpmtx[2, 2] = stsslmp[2]

        lmpmtx[1, 0] = stsslmp[3]
        lmpmtx[2, 0] = stsslmp[4]
        lmpmtx[2, 1] = stsslmp[5]

        lmpmtx[0, 1] = stsslmp[3]
        lmpmtx[0, 2] = stsslmp[4]
        lmpmtx[1, 2] = stsslmp[5]

        print basis * lmpmtx * basis.transpose()
        ###################################################################
        # vasp
        ###################################################################
        vaspcell = np.mat([[-1.7238,  1.5986,  1.5986],
                           [1.5986, -1.7238,  1.5986],
                           [1.7865,  1.7865, -1.5359]])
        invvasp = np.linalg.inv(vaspcell)
        invvasp /= np.linalg.norm(invvasp)
        stssvasp = np.array([-2.400844, -2.400844, 2.587547, -1.036506, 0.036871, 0.036871])

        vaspmtx[0, 0] = stssvasp[0]
        vaspmtx[1, 1] = stssvasp[1]
        vaspmtx[2, 2] = stssvasp[2]

        vaspmtx[1, 0] = stssvasp[3]
        vaspmtx[2, 0] = stssvasp[4]
        vaspmtx[2, 1] = stssvasp[5]

        vaspmtx[0, 1] = stssvasp[3]
        vaspmtx[0, 2] = stssvasp[4]
        vaspmtx[1, 2] = stssvasp[5]

        print basis * vaspmtx * basis.transpose()
        ###################################################################
        # vasp convention
        ###################################################################
        concell = np.mat([[3.4007,  0.0783,  0.0783],
                          [0.0783,  3.4007,  0.0783],
                          [-0.1566, -0.1566,  3.1658]])
        invcon = np.linalg.inv(concell)
        invcon /= np.linalg.norm(invcon)

        engy = -20.09869489
        stsscon = np.array([-2.51505, -2.51505, 3.934676, -1.279162, -0.002954, -0.002954])
        return

