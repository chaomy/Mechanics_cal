#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-26 00:16:08


import cal_md_gb_hcp
import cal_md_gb_pre
import cal_md_gb_pos
import cal_md_gb_run
import cal_md_gb_lmp
import cal_md_gb_indx
import cal_md_gb_hcp_0001
import cal_md_gb_hcp_1100
import cal_md_gb_hcp_1120
import get_data
import plt_drv
from md_pot_data import md_pot
from optparse import OptionParser
import gn_config


class md_gb(cal_md_gb_pre.md_gb_pre,
            cal_md_gb_run.md_gb_run,
            cal_md_gb_pos.md_gb_pos,
            cal_md_gb_lmp.md_gb_lmp,
            cal_md_gb_indx.md_gb_indx,
            cal_md_gb_hcp_0001.md_gb_loop,
            cal_md_gb_hcp.md_gb_hcp,
            cal_md_gb_hcp_1100.md_gb_hcp_1100,
            cal_md_gb_hcp_1120.md_gb_hcp_1120,
            plt_drv.plt_drv,
            get_data.get_data,
            gn_config.gnStructure):

    def __init__(self):
        # self.pot = md_pot.mg_kim
        self.pot = md_pot.mg_Poco
        cal_md_gb_pre.md_gb_pre.__init__(self)
        cal_md_gb_pos.md_gb_pos.__init__(self)
        cal_md_gb_run.md_gb_run.__init__(self)
        cal_md_gb_lmp.md_gb_lmp.__init__(self)
        cal_md_gb_indx.md_gb_indx.__init__(self)
        cal_md_gb_hcp_0001.md_gb_loop.__init__(self)
        cal_md_gb_hcp.md_gb_hcp.__init__(self)
        cal_md_gb_hcp_1100.md_gb_hcp_1100.__init__(self)
        cal_md_gb_hcp_1120.md_gb_hcp_1120.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = md_gb()
    dispatcher = {'gblist': drv.loop_gb_list,
                  'run': drv.loop_run,
                  'index': drv.hcp_tilt_index,
                  'build': drv.build_hcp_gb,
                  'plt': drv.loop_plt_angle,
                  '0001': drv.loop_angle_0001,
                  '1100': drv.loop_angle_1100,
                  '1120': drv.loop_angle_1120,
                  '1100g': drv.give_angle_1100}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
