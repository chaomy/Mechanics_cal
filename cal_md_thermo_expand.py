#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-27 01:26:23


from optparse import OptionParser
import matplotlib.pylab as plt
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import shutil
import gn_config
import get_data
import gn_lmp_infile
import gn_pbs
import plt_drv
import md_pot_data
import glob


class cal_md_thermo(gn_config.hcp,
                    gn_config.bcc,
                    gn_config.fcc,
                    get_data.get_data,
                    gn_pbs.gn_pbs,
                    gn_lmp_infile.gn_md_infile,
                    plt_drv.plt_drv):

    def __init__(self):
        self.indrv = gn_lmp_infile.gn_md_infile()
        self._pot = md_pot_data.md_pot.Nb_adp
        self._element = self._pot['element']
        self._pottype = self._pot['pair_style']
        self._size = np.array([16, 16, 16])
        self.set_lat('pot')
        self.unit_atoms = \
            ase.lattice.cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                            [0, 1, 0],
                                                            [0, 0, 1]],
                                                latticeconstant=self._lat,
                                                size=(1, 1, 1),
                                                symbol=self._element,
                                                pbc=(1, 1, 1))
        gn_config.bcc.__init__(self, self._pot)
        plt_drv.plt_drv.__init__(self)
        return

    def set_lat(self, tag='cal'):
        if tag is 'pot':
            self._lat = self._pot['lattice']
        elif tag is 'cal':
            os.system("lmp_mpi -i in.lat")
            lat = np.loadtxt("lat.txt")
            self._lat = np.average(lat)
        return

    def run_thermo(self, tag='run'):
        if tag == 'run':
            for i in range(1, 51):
                tend = 50 * i
                dirname = "dir-{:05.0f}".format(tend)
                os.chdir(dirname)
                os.system("mpirun -n 24 lmp_linux -i in.npt")
                os.chdir(os.pardir)
        elif tag == 'loopclc':
            temp_lx = []
            for i in range(1, 51):
                tend = 50 * i
                dirname = "dir-{:05.0f}".format(tend)
                os.chdir(dirname)
                temp_lx.append(self.get_temp_lat())
                os.chdir(os.pardir)
            print temp_lx
        temp_lx = np.array(temp_lx)
        np.savetxt("temp_lx.txt", temp_lx)
        return

    def draw_temp_vol(self):
        self.set_figs()
        datan = 10
        lastline = 31
        for i in range(1, 51):
            tend = 50 * i
            dirname = "dir-{:05.0f}".format(tend)
            print dirname
            print tend
            if os.path.isfile("{}/log.lammps".format(dirname)):
                raw = self.mreadlines("{}/log.lammps".format(dirname))
                for j in range(lastline, lastline + datan):
                    print raw[-j]
        return

    def given_temp_prep(self):
        unitatoms = self.unit_atoms.copy()
        atoms = unitatoms.repeat(([15, 15, 15]))
        self.write_lmp_config_data(atoms, "init.txt")
        oneatm = 1.01325
        tstart = 0.1
        initt = 50
        deltat = 50
        # initial
        dirname = "dir-{:05.0f}".format(0)
        self.mymkdir(dirname)
        self.indrv.write_md_thermo_expand('init')
        shutil.copy2("in.npt", dirname)
        # increase temp
        for i in range(50):
            # 2500  K
            tend = initt + deltat * i
            dirname = "dir-{:05.0f}".format(tend)
            self.mymkdir(dirname)
            if i == 0:
                self.indrv.write_md_thermo_expand(**{'tstart': tstart,
                                                     'tend': tend,
                                                     'pstart': 0.0,
                                                     'pend': oneatm})
            else:
                self.indrv.write_md_thermo_expand(**{'tstart': tstart,
                                                     'tend': tend,
                                                     'pstart': oneatm,
                                                     'pend': oneatm})
            tstart = tend
            shutil.copy2("in.npt", dirname)
        return

    def get_temp_lat(self, filename="log.lammps"):
        raw = self.mreadlines(filename)
        temp, lx, cnt = 0, 0, 0
        for i in range(len(raw)):
            line = raw[i]
            seg = line.split()
            if len(seg) > 0:
                if seg[0] == 'Step':
                    if int(raw[i + 1].split()[0]) > 1:
                        for j in range(522, 602):          # 600 lines
                            segs = raw[i + j].split()
                            temp += float(segs[3])
                            lx += float(segs[-1])
                            cnt += 1
                        break
        temp /= cnt
        lx /= cnt
        print temp, lx
        return (temp, lx)

    def theormo_expand_plt(self):
        temp_lx = np.loadtxt("temp_lx.txt")
        print temp_lx
        return

    def pressure_vs_vol(self, opt='prep'):
        lat0 = self._lat
        delta = -0.01
        npts = 30
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])
        lmp_bas = self.lmp_change_box(bas)

        if opt == 'clc':
            vol = np.zeros(npts)
            press = np.zeros(npts)

        for i in range(npts):
            rat = (1 + i * delta)**(1. / 3.)
            alat = rat * lat0
            dirname = 'dir-%04d' % (i)
            if opt == 'prep':
                self.mymkdir(dirname)
                cell = alat * lmp_bas
                atoms = ase.Atoms(self._pot['element'],
                                  positions=[[0, 0, 0]],
                                  cell=cell,
                                  pbc=[1, 1, 1])
                lmp_bas = self.lmp_change_box(bas)
                self.write_lmp_config_data(atoms, 'bulk.txt')
                os.system("mv bulk.txt  {}".format(dirname))
                os.system("cp in.init  {}".format(dirname))

            elif opt == 'run':
                os.chdir(dirname)
                os.system("lmp_mpi -i in.init")
                os.chdir(os.pardir)

            elif opt == 'clc':
                os.chdir(dirname)
                data = np.loadtxt("out.txt")
                print data
                vol[i] = data[0]
                press[i] = data[1]
                os.chdir(os.pardir)
                shutil.rmtree(dirname)

        if opt == 'clc':
            np.savetxt("data.txt", (vol, press))
        return

    def loop_pressure_vs_vol(self):
        dirlist = glob.glob("dir-*")
        if not os.path.isfile("in.lat"):
            os.system("cp ~/My_cal/Mechnical_cal/Bcc_Press_to_Vol/in.lat .")
            os.system("cp ~/My_cal/Mechnical_cal/Bcc_Press_to_Vol/in.init .")

        cnt = 0
        for mdir in dirlist[:]:
            print mdir
            if ((cnt % 1) == 0):
                if not os.path.isfile("fig-{}.png".format(mdir)):
                    os.system("cp {}/dummy.lammps.ADP .".format(mdir))
                    self.set_lat('cal')
                    self.pressure_vs_vol('prep')
                    self.pressure_vs_vol('run')
                    self.pressure_vs_vol('clc')
                    self.vasp_energy_stress_vol_plt(mdir, 30)
                    os.remove("dummy.lammps.ADP")
                    os.system("mv data.txt  {}".format(mdir))
                    os.system("cp p2v.png fig-{}.png".format(mdir))
            cnt += 1
        return

    def delete_some_cands(self):
        raw = self.mreadlines("candlist.txt")
        self.mymkdir("save")
        for line in raw:
            figname = "fig-{}.png".format(line.split()[0])
            os.system("mv {} save".format(figname))
        return

    def vasp_energy_stress_vol_plt(self,
                                   inlabel='p-v',
                                   npt=30):
        self.set_keys("upper right")
        self.set_111plt((10, 6.5))
        (vol, press) = np.loadtxt("data.txt")
        (dft_vol, dft_press) = np.loadtxt("../../DATA_DFT_PV.txt")

        vol = vol / vol[0]
        vol = vol**3
        dft_vol = dft_vol / dft_vol[0]

        self.ax.plot(vol[:npt], press[:npt],
                     label='adp',
                     **next(self.keysiter))

        self.ax.plot(dft_vol[:npt], dft_press[:npt],
                     label='pbe',
                     **next(self.keysiter))
        plt.xlabel('relative volume (V / V$_0$)',
                   {'fontsize': self.myfontsize})
        plt.ylabel('presssure (GPa)', {'fontsize': self.myfontsize})
        self.add_legends(self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig("p2v.png", **self.figsave)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_md_thermo()
    if options.mtype.lower() == 'thermo':
        drv.given_temp_prep()

    if options.mtype.lower() == 'runjob':
        drv.run_thermo()

    if options.mtype.lower() == 'clc':
        #  drv.get_temp_lat()
        drv.run_thermo('loopclc')

    if options.mtype.lower() == 'p2vprep':
        drv.pressure_vs_vol('prep')

    if options.mtype.lower() == 'p2vrun':
        drv.pressure_vs_vol('run')

    if options.mtype.lower() == 'p2vclc':
        drv.pressure_vs_vol('clc')

    if options.mtype.lower() == 'p2vauto':
        if not os.path.isfile("in.init"):
            os.system(
                "cp ~/My_cal/Mechnical_cal/Bcc_Press_to_Vol/in.lat_adp in.lat")
            os.system(
                "cp ~/My_cal/Mechnical_cal/Bcc_Press_to_Vol/in.init_adp in.init")
        drv.set_lat('cal')
        drv.pressure_vs_vol('prep')
        drv.pressure_vs_vol('run')
        drv.pressure_vs_vol('clc')
        drv.vasp_energy_stress_vol_plt(30)

    if options.mtype.lower() == 'p2vplt':
        drv.vasp_energy_stress_vol_plt(30)

    if options.mtype.lower() == 'plot':
        drv.draw_temp_vol()

    if options.mtype.lower() == 'expandplt':
        drv.theormo_expand_plt()

    if options.mtype.lower() == 'loopp2v':
        drv.loop_pressure_vs_vol()

    if options.mtype.lower() == 'del':
        drv.delete_some_cands()
