#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-07 00:49:31


from copy import deepcopy
import os
import numpy as np
import md_pot_data
import ase.io
import gsf_data


class cal_qe_gsf(object):

    def __init__(self):
        self.vcaWTa = {'WTa0.05': md_pot_data.qe_pot.vca_W95Ta05,
                       'WTa0.10': md_pot_data.qe_pot.vca_W90Ta10,
                       'WTa0.15': md_pot_data.qe_pot.vca_W85Ta15,
                       'WTa0.20': md_pot_data.qe_pot.vca_W80Ta20,
                       'WTa0.25': md_pot_data.qe_pot.vca_W75Ta25}
        self.vcapots = {
            'WTa50': md_pot_data.qe_pot.vca_W50Ta50,
            'WRe00': md_pot_data.qe_pot.pbe_w,
            'WRe05': md_pot_data.qe_pot.vca_W95Re05,
            'WRe10': md_pot_data.qe_pot.vca_W90Re10,
            'WRe15': md_pot_data.qe_pot.vca_W85Re15,
            'WRe20': md_pot_data.qe_pot.vca_W80Re20,
            'WRe25': md_pot_data.qe_pot.vca_W75Re25,
            'WRe50': md_pot_data.qe_pot.vca_W50Re50}
        return

    def setup_qe_scf(self):
        self.set_cal_type('scf')
        self.set_ecut('41')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-6')
        self.set_kpnts(gsf_data.gsfkpts[self.mgsf])
        self.set_maxseconds(3600 * 24)
        return

    def setup_qe_relax(self):
        self.set_cal_type('relax')
        self.set_ecut('40')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-6')
        self.set_kpnts(gsf_data.gsfkpts[self.mgsf])
        self.set_maxseconds(3600 * 100)
        return

    def gn_qe_single_dir_gsf(self, mtype='relax'):
        if self.mgsf in ['x111z110']:
            atoms = self.gn_bcc110()
        elif self.mgsf in ['x111z112']:
            atoms = self.gn_gsf_atoms()
        atoms.wrap()
        perf_cells = deepcopy(atoms.get_cell())
        ase.io.write('perf_poscar', images=atoms, format='vasp')

        # disps = np.arange(0.42, 0.58, 0.04)
        disps = np.arange(0.34, 0.48, 0.04)
        disps = np.append(disps, 0.0)
        npts = len(disps)

        if mtype in ['scf']:
            self.setup_qe_scf()
        elif mtype in ['relax']:
            self.setup_qe_relax()

        for i, disp in zip(range(npts), disps):
            dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            self.mymkdir(dirname)
            os.chdir(dirname)
            if self.mgsf in ['x111z110']:
                disp_vector = [disp, disp, 0]
            else:
                disp_vector = [disp, 0.0, 0.0]
            disp_matrix_direct = self.gn_displacement(
                atoms.copy(), disp_vector)
            disp_matrix = deepcopy(disp_matrix_direct)
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * perf_cells[0, 0]
            print disp_matrix
            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            if mtype in ['scf']:
                self.gn_qe_scf_tf(local_atoms)
            elif mtype in ['relax']:
                self.gn_qe_relax_tf(local_atoms, 'xy')

            self.set_pbs(dirname)
            os.system("cp $POTDIR/{} . ".format(self.pot['file']))
            ase.io.write('poscar', images=local_atoms, format='vasp')
            self.write_lmp_config_data(local_atoms)
            os.system("mv poscar ../poscar.{:03d}".format(i))
            os.chdir(os.pardir)
        return

    def check_gsf(self):
        for key in gsf_data.gsfsize:
            atomss = self.set_bcc_convention(
                in_direction=gsf_data.gsfbase[key],
                in_size=gsf_data.gsfsize[key])

            atomsb = self.set_bcc_convention(
                in_direction=gsf_data.gsfbase[key],
                in_size=gsf_data.bulksize[key])

            print(key, len(atomss) - len(atomsb),
                  gsf_data.gsfpopn[key])
        return
