#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-28 22:31:58


from optparse import OptionParser
from collections import OrderedDict
import os
import gsf_data
import cal_qe_gsf
import numpy as np
import md_pot_data
import ase.io


class cal_surface(cal_qe_gsf.cal_gsf):

    def __init__(self,
                 pot=md_pot_data.qe_pot.pbe_w,
                 msurf='x100z100'):
        self.pot = pot
        self.msurf = msurf
        cal_qe_gsf.cal_gsf.__init__(self, self.pot, self.msurf)
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)
        return

    def gn_surface_atoms(self):
        mgsf = self.mgsf
        atomss = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
            in_size=gsf_data.gsfsize[mgsf])
        for i in range(gsf_data.gsfpopn[mgsf]):
            atomss.pop()

        atomsb = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
            in_size=gsf_data.bulksize[mgsf])
        atomss.wrap()
        atomsb.wrap()
        return [atomss, atomsb]

    def prep_qe_surface(self, dirtag='dir'):
        configs = ['surf', 'bulk']
        atomsl = self.gn_surface_atoms()
        self.setup_qe_scf()
        for config, atoms in zip(configs, atomsl):
            dirname = '{}-{}-{}'.format(
                dirtag, self.mgsf, config)
            self.mymkdir(dirname)

            os.chdir(dirname)
            self.gn_qe_scf_tf(atoms)
            ase.io.write('poscar.vasp', images=atoms, format='vasp')
            os.system("cp $POTDIR/{} . ".format(self.pot['file']))
            self.set_pbs(dirname)
            os.chdir(os.pardir)
        return

    def cal_qe_surface(self, dirtag='dir'):
        configs = ['bulk', 'surf']
        data = np.zeros(len(configs))
        for i, config in zip(range(len(configs)), configs):
            dirname = '{}-{}-{}'.format(
                dirtag, self.mgsf, config)
            os.chdir(dirname)
            self.qe_get_cell()
            area = self.cal_xy_area()
            data[i] = self.qe_get_energy_stress()[0]
            os.chdir(os.pardir)
        return np.array([area, data[0], data[1]])

    def loop_pot_surf(self, intag='prep'):
        vcapots = {
            'WRe00': md_pot_data.qe_pot.pbe_w,
            'WRe05': md_pot_data.qe_pot.vca_W95Re05,
            'WRe10': md_pot_data.qe_pot.vca_W90Re10,
            'WRe15': md_pot_data.qe_pot.vca_W85Re15,
            'WRe20': md_pot_data.qe_pot.vca_W80Re20,
            'WRe25': md_pot_data.qe_pot.vca_W75Re25,
            'WRe50': md_pot_data.qe_pot.vca_W50Re50}
        surfs = ['x100z100']
        for key in vcapots:
            for surf in surfs:
                dirtag = 'dir-{}'.format(key)
                self.__init__(vcapots[key], surf)
                self.prep_qe_surface(dirtag)
        return

    def loop_cal_surf(self):
        vcapots = OrderedDict([
            ('WRe00', md_pot_data.qe_pot.pbe_w),
            ('WRe05', md_pot_data.qe_pot.vca_W95Re05),
            ('WRe10', md_pot_data.qe_pot.vca_W90Re10),
            ('WRe15', md_pot_data.qe_pot.vca_W85Re15),
            ('WRe20', md_pot_data.qe_pot.vca_W80Re20),
            ('WRe25', md_pot_data.qe_pot.vca_W75Re25),
            ('WRe50', md_pot_data.qe_pot.vca_W50Re50)])
        surf = 'x100z100'
        npts = len(vcapots)
        data = np.ndarray([npts, 5])
        for key, i in zip(vcapots.keys(), range(npts)):
            dirtag = 'dir-{}'.format(key)
            self.__init__(vcapots[key], surf)
            if i == 6:
                data[i, 0] = 0.5
            else:
                data[i, 0] = 0.05 * i
            data[i, 2:] = self.cal_qe_surface(dirtag)
            data[i, 1] = 0.5 * ((data[i, -1] - data[i, -2]) / data[i, -3])
        np.savetxt('surf.dat', data)
        return

    def loop_plt_surf(self):
        np.loadtxt()
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    (options, args) = parser.parse_args()
    drv = cal_surface()
    dispatcher = {'prep': drv.prep_qe_surface,
                  'chk': drv.check_gsf,
                  'loopprep': drv.loop_pot_surf,
                  'loopsurf': drv.loop_cal_surf}
    dispatcher[options.mtype.lower()]()
