#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-08 21:10:27


from optparse import OptionParser
from scipy.optimize import leastsq
import plt_drv
import glob
import os
import numpy as np
import md_pot_data
import gn_config
import gn_qe_inputs
import get_data
import gn_kpoints
import gn_incar
import gn_pbs
import output_data


class cal_cij(gn_config.bcc,
              gn_config.fcc,
              gn_config.hcp,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              output_data.output_data,
              gn_qe_inputs.gn_qe_infile,
              plt_drv.plt_drv):

    def __init__(self, inpot=md_pot_data.qe_pot.pbe_w):
        self.unit_delta = 0.001
        self.looptimes = 5 
        self.volume = None
        self.energy0 = None
        self.pot = inpot
        self.cij_type_list = ['c11', 'c12', 'c44']
        self.cij_type = 'c11'

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        output_data.output_data.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)

        if self.pot['structure'] == 'bcc':
            gn_config.bcc.__init__(self, self.pot)
        elif self.pot['structure'] == 'fcc':
            gn_config.fcc.__init__(self, self.pot)
        elif self.pot['structure'] == 'hcp':
            gn_config.hcp.__init__(self, self.pot)
        return

    def fit_para(self, delta_list, energy_list):
        xdata = delta_list
        ydata = energy_list - self.energy0

        def residuals(params):
            a, c = params
            return ydata - (0.5 * a * (xdata**2) + c)
        r = leastsq(residuals, [1.0, 0.0])
        print r
        a, c = r[0]
        print "vol = ", self.volume, "A^3"
        a = a / self.volume * self.ev_angstrom3_to_GPa 
        return a

    def obtain_cij_old(self):
        self.volume = 3.071**3
        del2coeffs = np.transpose(np.zeros(3))
        fileList = ['c11_summary', 'c12_summary', 'c44_summary']
        for file, i in zip(fileList, range(3)):
            raw = np.loadtxt(file)
            self.energy0 = raw[0, 1] * 13.605698066
            delta_list, energy_list = raw[:, 0], raw[:, 1] * 13.605698066
            del2coeffs[i] = self.fit_para(delta_list, energy_list)
        # convmat = np.mat([[3, 6, 0], [2, -2, 0], [0, 0, 4]])
        # print 2 * np.linalg.pinv(convmat) * np.transpose(np.mat(del2coeffs))
        return

    def obtain_cij(self, opt='np'):
        del2coeffs = np.transpose(np.zeros(3))
        self.set_volume_energy0()
        filelist = ['data_c11.txt', 'data_c12.txt', 'data_c44.txt']
        for i in range(3):
            raw = np.loadtxt(filelist[i])
            delta_list, energy_list = raw[:, 0], raw[:, 1]
            del2coeffs[i] = self.fit_para(delta_list, energy_list)

        # print np.linalg.pinv(convmat) * np.transpose(np.mat(del2coeffs))
        print del2coeffs
        c11 = (0.111111111111111 * del2coeffs[0] + 0.333333333333333 * del2coeffs[1])
        c12 = (0.111111111111111 * del2coeffs[0] - 0.166666666666667 * del2coeffs[1])
        c44 = 0.25 * del2coeffs[2]
        print(c11, c12, c44)
        # with open("cij.dat", 'w') as fout:
        #     fout.write("C11\t%f\t\nC12\t%f\t\nC44\t%f\t\n" % (c11, c12, c44))
        #     fout.close()
        return

    def set_volume_energy0(self):
        (engy, vol, stress) = self.qe_get_energy_stress(
            filename='dir-c11-p000/qe.out')
        self.volume = vol * md_pot_data.unitconv.ulength["BohrtoA"]**3
        print self.volume
        self.energy0 = engy
        return

    def set_cij_type(self, cij_type):
        self.cij_type = cij_type
        return

    def gn_qe_cij_infile(self, atoms):
        self.set_thr('1.0D-6')
        self.set_ecut('40')
        self.set_kpnts((39, 39, 39))
        # self.set_kpnts((19, 19, 19))
        self.set_degauss('0.03D0')
        with open('qe.in', 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid)
            fid.close()
        return

    def loop_prepare_cij(self):
        for i in range(len(self.cij_type_list)):
            mtype = self.cij_type_list[i]
            self.set_cij_type(mtype)
            for j in range(-self.looptimes, self.looptimes):
                delta = self.unit_delta * j
                if j >= 0:
                    dirname = "dir-%s-p%03d" % (mtype, j)
                else:
                    dirname = "dir-%s-n%03d" % (mtype, np.abs(j))
                self.mymkdir(dirname)
                os.chdir(dirname)
                self.set_pbs(dirname)
                # atoms = self.write_bcc_primitive_with_strain(delta=delta,
                #                                              in_tag=mtype,
                #                                              write=False)
                atoms = self.write_bcc_with_strain(delta=delta,
                                                   in_tag=mtype,
                                                   write=False)
                self.gn_qe_cij_infile(atoms)
                os.system("cp $POTDIR/{} .".format(self.pot['file']))
                os.chdir(os.pardir)
        return

    def cij_plt(self):
        filelist = ['data_c11.txt', 'data_c12.txt', 'data_c44.txt']
        for file in filelist:
            raw = np.loadtxt(file)
            delta_list, energy_list = raw[:, 0], raw[:, 1]
            self.set_111plt()
            self.ax.plot(delta_list, energy_list,
                         label='cij', **next(self.keysiter))
            self.fig.savefig('fig_{}.png'.format(file[:-8:-4]), **self.figsave)
        return

    def collect_data_cij(self):
        for i in range(len(self.cij_type_list)):
            mtype = self.cij_type_list[i]
            self.set_cij_type(mtype)
            out_file_name = "data_%s.txt" % (mtype)
            for j in range(-self.looptimes, self.looptimes):
                delta = self.unit_delta * j
                if j >= 0:
                    dirname = "dir-%s-p%03d" % (mtype, j)
                else:
                    dirname = "dir-%s-n%03d" % (mtype, np.abs(j))
                os.chdir(dirname)
                print "i am in ", dirname
                (energy, vol, stress) = self.qe_get_energy_stress()  # in eV
                os.chdir(os.pardir)
                self.output_delta_energy(delta,
                                         energy,
                                         file_name=out_file_name)
        return

    def set_pbs(self, dirname):
        self.set_wall_time(2)
        self.set_job_title(dirname)
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun pw.x < qe.in > qe.out")
        self.write_pbs(od=True)
        return

    def loop_sub_drvs(self):
        dir_list = glob.glob("dir-*")
        for i in range(len(dir_list)):
            os.chdir(dir_list[i])
            os.system("qsub va.pbs")
            os.chdir(os.pardir)
        return

    def loop_pots(self):
        potlist = {'WTa0.25': md_pot_data.qe_pot.vca_W75Ta25,
                   'WTa0.20': md_pot_data.qe_pot.vca_W80Ta20,
                   'WTa0.15': md_pot_data.qe_pot.vca_W85Ta15,
                   'WTa0.10': md_pot_data.qe_pot.vca_W90Ta10,
                   'WTa0.05': md_pot_data.qe_pot.vca_W95Ta05}
        for key in potlist.keys():
            self.mymkdir(key)
            self.__init__(potlist[key])
            os.chdir(key)
            self.loop_prepare_cij()
            os.chdir(os.pardir)
        return


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_cij()
    dispatcher = {'pots': drv.loop_pots,
                  'prep': drv.loop_prepare_cij,
                  'sub': drv.loop_sub_drvs,
                  'cij': drv.obtain_cij,  # obtain_cij_old
                  'clc': drv.collect_data_cij,
                  'plt': drv.cij_plt}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

    #  A.set_volume_energy0()
    #  A.obtain_cij_old()
    #  c12 = 0.25 * C11plus2C12 + 0.75 * C11minusC12
    #  c11 = 0.25 * C11plus2C12 - 0.25 * C11minusC12
    #  c44 = float(get_cxx("c44",Lattice_constant,delta_increment,nt,E0)) * (1./4.) * inv_V0  * unit_conversion1
