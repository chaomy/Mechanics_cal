#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-09 01:51:07


from optparse import OptionParser
from collections import OrderedDict
from md_pot_data import fluxdirs
import os
import numpy as np
import md_pot_data
import ase.io

vcapotsRe = OrderedDict([
    ('WRe00', md_pot_data.qe_pot.pbe_w),
    ('WRe05', md_pot_data.qe_pot.vca_W95Re05),
    ('WRe10', md_pot_data.qe_pot.vca_W90Re10),
    ('WRe15', md_pot_data.qe_pot.vca_W85Re15),
    ('WRe20', md_pot_data.qe_pot.vca_W80Re20),
    ('WRe25', md_pot_data.qe_pot.vca_W75Re25),
    ('WRe50', md_pot_data.qe_pot.vca_W50Re50),
    ('WTa50', md_pot_data.qe_pot.vca_W50Ta50)])

vcapots = OrderedDict([
    ('WTa05', md_pot_data.qe_pot.vca_W95Ta05),
    ('WTa10', md_pot_data.qe_pot.vca_W90Ta10),
    ('WTa15', md_pot_data.qe_pot.vca_W85Ta15),
    ('WTa20', md_pot_data.qe_pot.vca_W80Ta20),
    ('WTa25', md_pot_data.qe_pot.vca_W75Ta25),
    ('WTa50', md_pot_data.qe_pot.vca_W50Ta50)])

class cal_qe_surface(object):

    def prep_qe_surface(self, dirtag='dir', mtype='relax'):
        configs = ['surf', 'bulk']
        atomsl = self.gn_surface_atoms()

        if mtype in ['scf']:
            self.setup_qe_scf()
        elif mtype in ['relax']:
            self.setup_qe_relax()

        for config, atoms in zip(configs, atomsl):
            mdir = '{}-{}-{}'.format(dirtag, self.mgsf, config)
            self.mymkdir(mdir)

            os.chdir(mdir)

            if mtype in ['relax']:
                self.gn_qe_relax_tf(atoms)
            elif mtype in ['scf']:
                self.gn_qe_scf_tf(atoms)

            ase.io.write('poscar.vasp', images=atoms, format='vasp')
            os.system("cp $POTDIR/{} . ".format(self.pot['file']))
            self.set_pbs(mdir, 'qe')
            os.chdir(os.pardir)
        return

    def cal_qe_surface(self, dirtag='dir'):
        configs = ['bulk', 'surf']
        data = np.zeros(len(configs))
        for i, config in zip(range(len(configs)), configs):
            mdir = '{}-{}-{}'.format(
                dirtag, self.mgsf, config)
            os.chdir(mdir)
            self.qe_get_cell()
            area = self.cal_xy_area()
            data[i] = self.qe_get_energy_stress()[0]
            os.chdir(os.pardir)
        return np.array([area, data[0], data[1]])

    def loop_pot_surf(self, intag='prep'):
        surfs = ['x100z100']
        for key in vcapots:
            for surf in surfs:
                dirtag = 'dir-{}'.format(key)
                self.__init__(vcapots[key], surf)
                self.prep_qe_surface(dirtag)
        return

    def loop_cal_surf(self):
        surf = 'x100z100'
        npts = len(vcapots)
        data = np.ndarray([npts, 5])
        ev2j = 1.60218
        for key, i in zip(vcapots.keys(), range(npts)):
            dirtag = 'dir-{}'.format(key)
            self.__init__(vcapots[key], surf)
            if i == 5:
                data[i, 0] = 0.5
            else:
                data[i, 0] = 0.05 * i
            data[i, 2:] = self.cal_qe_surface(dirtag)
            data[i, 1] = 0.5 * \
                ((data[i, -1] - data[i, -2]) / data[i, -3]) * ev2j
        np.savetxt('surf.dat', data)
        return

    def plt_surf(self):
        data = np.loadtxt('surf.dat')
        self.set_111plt()
        axlist = [self.ax]
        self.ax.plot(data[:, 0], data[:, 1],
                     label='(100)', **next(self.keysiter))
        self.add_legends(*axlist)
        self.fig.savefig('fig_surfE.png')
        return

    def transdata(self, ptype='scp'):
        npts = len(vcapots)
        configs = ['surf', 'bulk']
        for config in configs:
            for key, i in zip(vcapots.keys(), range(npts)):
                dirtag = 'dir-{}'.format(key)
                mdir = '{}-{}-{}'.format(
                    dirtag, self.mgsf, config)
                print mdir
                self.mymkdir(mdir)
                fdir = fluxdirs['QE'] + \
                    'VC_WTa/surf/{}'.format(mdir)
                os.system('scp {}/qe.out {}'.format(fdir, mdir, mdir))
                os.system('scp {}/qe.in {}'.format(fdir, mdir, mdir))

        # fdir = fluxdirs['QE'] + \
        #     'VC_WRe/Bcc_QE_VCA_WRe_unrelaxgsf/{}/'.format(
        #         dirtree[self.mgsf][tag])
        # os.system('scp {}/{}/qe.out {}'.format(fdir, mdir, mdir))
        # os.system('scp {}/{}/qe.in {}'.format(fdir, mdir, mdir))
        return
