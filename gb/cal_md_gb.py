#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-21 09:44:17


from get_data import get_data
from cal_md_gb_hcp import md_gb_hcp
from cal_md_gb_pre import md_gb_pre
from cal_md_gb_pos import md_gb_pos
from cal_md_gb_run import md_gb_run
from cal_md_gb_lmp import md_gb_lmp
from cal_md_gb_indx import md_gb_indx
from cal_md_gb_hcp_0001 import md_gb_loop
from cal_md_gb_hcp_1100 import md_gb_hcp_1100 
from cal_md_gb_hcp_1120 import md_gb_hcp_1120 
from plt_drv import plt_drv
from md_pot_data import md_pot
from optparse import OptionParser
from gn_config import gnStructure


class md_gb(md_gb_pre,
            md_gb_run,
            md_gb_pos,
            md_gb_lmp,
            md_gb_indx,
            md_gb_loop,
            md_gb_hcp,
            md_gb_hcp_1100,
            md_gb_hcp_1120,
            plt_drv,
            get_data,
            gnStructure):

    def __init__(self):
        # self.pot = md_pot.mg_kim
        self.pot = md_pot.mg_Poco
        md_gb_pre.__init__(self)
        md_gb_pos.__init__(self)
        md_gb_run.__init__(self)
        md_gb_lmp.__init__(self)
        md_gb_indx.__init__(self)
        md_gb_loop.__init__(self)
        md_gb_hcp.__init__(self)
        md_gb_hcp_1100.__init__(self)
        md_gb_hcp_1120.__init__(self)
        plt_drv.__init__(self)
        get_data.__init__(self)
        gnStructure.__init__(self, self.pot)


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
                  'angle': drv.loop_angle,
                  'plt': drv.loop_plt_angle,
                  '1100': drv.loop_angle_1100,
                  '1120': drv.loop_angle_1120}

    # 'loop': drv.loop_dispx_hcp,
    # 'thk': drv.loop_thickness

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
