#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-20 13:27:31


from optparse import OptionParser
import ase
import ase.io
import copy
import os
import numpy as np
import glob

try:
    import atomman as am
    import atomman.lammps as lmp
    import atomman.unitconvert as uc
    import gn_config
    import get_data
    import gn_lmp_infile
    import gn_pbs

except ImportError:
    print("error during import")


class cal_md_hcp_basic(gn_config.hcp,
                       get_data.get_data,
                       gn_pbs.gn_pbs,
                       gn_lmp_infile.gn_md_infile):

    def __init__(self):
        gn_lmp_infile.gn_md_infile.__init__(self)
        return

    def hcp_lattice(self):
        self.write_hcp_lattice_infile("in.hcp",
                                      2.85)
        return

if __name__ == '__main__':
    drv = cal_md_hcp_basic()
    drv.hcp_lattice()
