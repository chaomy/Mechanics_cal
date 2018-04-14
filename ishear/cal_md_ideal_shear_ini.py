#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-05 02:31:51


from cal_md_ideal_shear_pre import cal_bcc_ideal_shear_pre
from cal_md_ideal_shear_run import cal_bcc_ideal_shear_run
from cal_md_ideal_shear_pos import cal_bcc_ideal_shear_pos
from cal_md_ideal_shear_plt import cal_bcc_ideal_shear_plt
from cal_md_ideal_shear_mesh import cal_bcc_ideal_shear_mesh
import os
import numpy as np
import gn_config
import gn_pbs
import get_data
import plt_drv
import gn_qe_inputs


class cal_bcc_ideal_shear(get_data.get_data,
                          gn_config.gnStructure,
                          gn_pbs.gn_pbs,
                          plt_drv.plt_drv,
                          gn_qe_inputs.gn_qe_infile,
                          cal_bcc_ideal_shear_pre,
                          cal_bcc_ideal_shear_run,
                          cal_bcc_ideal_shear_pos,
                          cal_bcc_ideal_shear_plt,
                          cal_bcc_ideal_shear_mesh):

    def __init__(self, inpot, shtype='110'):
        self.pot = self.load_data('../BASICS/pot.dat')
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        cal_bcc_ideal_shear_pre.__init__(self)
        cal_bcc_ideal_shear_run.__init__(self)
        cal_bcc_ideal_shear_pos.__init__(self)
        cal_bcc_ideal_shear_plt.__init__(self)
        cal_bcc_ideal_shear_mesh.__init__(self)
        self.npts = 20
        self.delta = 0.02
        shd111p211 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., 1., -2.]),
                      'e3': np.array([-1., 1., 0.])}
        shd111p110 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., -1, 0]),
                      'e3': np.array([1, 1., -2])}
        if shtype == '211':
            e1 = shd111p211['e1']
            e2 = shd111p211['e2']
            e3 = shd111p211['e3']
        elif shtype == '110':
            e1 = shd111p110['e1']
            e2 = shd111p110['e2']
            e3 = shd111p110['e3']
        e1 = e1 / np.linalg.norm(e1)
        e2 = e2 / np.linalg.norm(e2)
        e3 = e3 / np.linalg.norm(e3)
        self.basis = np.mat([e1, e2, e3])
        self.setup_qe_params()