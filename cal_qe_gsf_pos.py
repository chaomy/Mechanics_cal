#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-09 21:07:12


from itertools import cycle
from numpy import arange, min, max, append
from numpy import loadtxt, savetxt, ndarray
from md_pot_data import fluxdirs
import os


# dirtree = {'x111z110': {
#     '00': 'Bcc_QE_VCA_WRe00_gsfx111z110',
#     '05': 'Bcc_QE_VCA_WRe05_gsfx111z110',
#     '10': 'Bcc_QE_VCA_WRe10_gsfx111z110',
#     '15': 'Bcc_QE_VCA_WRe15_gsfx111z110',
#     '20': 'Bcc_QE_VCA_WRe20_gsfx111z110',
#     '25': 'Bcc_QE_VCA_WRe25_gsfx111z110',
#     '50': 'Bcc_QE_VCA_WRe50_gsfx111z110'
# }, 'x111z112': {
#     '00': 'Bcc_QE_VCA_WRe00_gsfx111z112',
#     '05': 'Bcc_QE_VCA_WRe05_gsfx111z112',
#     '10': 'Bcc_QE_VCA_WRe10_gsfx111z112',
#     '15': 'Bcc_QE_VCA_WRe15_gsfx111z112',
#     '20': 'Bcc_QE_VCA_WRe20_gsfx111z112',
#     '25': 'Bcc_QE_VCA_WRe25_gsfx111z112',
#     '50': 'Bcc_QE_VCA_WRe50_gsfx111z112'}
# }

dirtree = {'x111z110': {
    '00': 'Bcc_WRe00_gsfx111z110',
    '05': 'Bcc_WRe05_gsfx111z110',
    '10': 'Bcc_WRe10_gsfx111z110',
    '15': 'Bcc_WRe15_gsfx111z110',
    '20': 'Bcc_WRe20_gsfx111z110',
    '25': 'Bcc_WRe25_gsfx111z110',
    '50': 'Bcc_WRe50_gsfx111z110',
    'ta': 'Bcc_WTa50_gsfx111z110'
}, 'x111z112': {
    '00': 'Bcc_WRe00_gsfx111z112',
    '05': 'Bcc_WRe05_gsfx111z112',
    '10': 'Bcc_WRe10_gsfx111z112',
    '15': 'Bcc_WRe15_gsfx111z112',
    '20': 'Bcc_WRe20_gsfx111z112',
    '25': 'Bcc_WRe25_gsfx111z112',
    '50': 'Bcc_WRe50_gsfx111z112',
    'ta': 'Bcc_WTa50_gsfx111z112'}
}


class cal_qe_gsf_pos(object):

    def loop_clcenergy(self):
        mlist = ['00', '05', '10', '15', '20', '25', '50']
        for mkey in mlist:
            mdir = dirtree[self.mgsf][mkey]
            os.chdir(mdir)
            self.clc_qe_gsf_engy()
            os.system('mv gsf.dat ../dat.{}.txt'.format(mdir))
            os.chdir(os.pardir)
        return

    def transdata(self, ptype='scp', tag='ta'):
        # disps = arange(0.46, 0.56, 0.04)
        disps = arange(0.42, 0.48, 0.04)
        disps = append(disps, 0.0)
        for disp in disps:
            mdir = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            self.mymkdir(mdir)
            if ptype in ['scp']:
                fdir = fluxdirs['QE'] + \
                    'VC_WRe/Bcc_QE_VCA_WRe_relaxgsf/{}'.format(
                        dirtree[self.mgsf][tag])
                os.system('scp {}/{}/qe.out {}'.format(fdir, mdir, mdir))
                os.system('scp {}/{}/qe.in {}'.format(fdir, mdir, mdir))
                print fdir
        return

    def clc_qe_gsf_engy(self, fname='gsf'):
        disps = 0.0
        # disps = append(disps, arange(0.46, 0.56, 0.04))
        disps = append(disps, arange(0.42, 0.48, 0.04))
        # disps = append(disps, 1.0)
        npts = len(disps)
        data = ndarray([npts, 4])
        for i, disp in zip(range(npts), disps):
            if disp == 1:
                dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, 0.0)
            else:
                dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            os.chdir(dirname)
            print dirname
            # print(self.qe_get_cell())
            data[i, 0] = i
            data[i, 1] = disp
            data[i, 2] = self.cal_xy_area()
            data[i, 3] = self.qe_get_energy_stress()[0]
            os.chdir(os.pardir)
        savetxt('gsf.dat', data)
        return

    def loop_plt_gsf(self):
        mlist = ['00', '05', '10', '15', '20', '25', '50']
        lablist = ['W', r'WRe$_{0.05}$', r'WRe$_{0.10}$',
                   r'WRe$_{0.15}$', r'WRe$_{0.20}$',
                   r'WRe$_{0.25}$', r'WRe$_{0.50}$']
        toldat = ndarray([len(mlist), 2])
        self.set_111plt()
        self.set_keys(ncol=3, lg='in')
        axlist = [self.ax]
        reratio = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50]
        for mkey, labtag, re, i in zip(mlist, lablist,
                                       reratio, range(len(reratio)))[:6]:
            mdir = dirtree[self.mgsf][mkey]
            mfile = 'dat.{}.txt'.format(mdir)
            data = loadtxt(mfile)
            gsf = (data[:, 3] - min(data[:, 3])) / (2. * data[:, 2])
            self.ax.plot(data[:, 1], gsf,
                         label=labtag, **next(self.keysiter))
            toldat[i, 0] = re
            toldat[i, 1] = max(gsf)
        ylabiter = cycle([r'Stacking fault energy [$eV/A^2$]'])
        xlabiter = cycle(['Normalized displacement along [111](211)'])
        self.ax.set_ylim([-0.2e-2, 8.15e-2])
        self.add_y_labels(ylabiter, *axlist)
        self.add_x_labels(xlabiter, *axlist)
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.fig.savefig('gsf_{}.png'.format(self.mgsf), **self.figsave)
        savetxt('tol_{}.txt'.format(self.mgsf), toldat)
        return

    def plt_tol(self):
        dat1 = loadtxt('tol_x111z110.txt')[:6]
        dat2 = loadtxt('tol_x111z112.txt')[:6]
        self.set_111plt()
        axlist = [self.ax]
        coeff = 16.021766208
        print coeff * dat1[:, 1]
        print coeff * dat2[:, 1]
        self.ax.plot(dat1[:, 0], dat1[:, 1],
                     label='USFE [111](110)', **next(self.keysiter))
        self.ax.plot(dat2[:, 0], dat2[:, 1],
                     label='USFE [111](112)', **next(self.keysiter))
        ylabiter = cycle([r'Stacking fault energy[$eV/A^2$]'])
        xlabiter = cycle(['Re concentration'])
        self.add_y_labels(ylabiter, *axlist)
        self.add_x_labels(xlabiter, *axlist)
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.fig.savefig('fig_tol.png'.format(self.mgsf), **self.figsave)
        self.closefig()

        # plot ductility parameters
        surfE03 = loadtxt('surfE03.txt')
        print 'SurfE'
        print surfE03 * coeff
        self.set_111plt()
        self.ax.plot(dat1[:, 0], surfE03 / dat1[:, 1],
                     label='D [111](110)', **next(self.keysiter))
        self.ax.plot(dat1[:, 0], surfE03 / dat2[:, 1],
                     label='D [111](112)', **next(self.keysiter))
        self.add_legends(self.ax)
        self.add_x_labels(cycle(['Re concentration']), self.ax)
        self.add_y_labels(cycle(['Ductility parameter']), self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig('fig_duc.png', **self.figsave)
        return
