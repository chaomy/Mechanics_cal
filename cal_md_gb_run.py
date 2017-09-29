#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-29 02:11:54

import os
from glob import glob
from subprocess import Popen

class md_gb_run(object):

    def loop_run(self):
        dlist = glob("gb-*")
        for mdir in dlist:
            os.chdir(mdir)
            sublist = glob("mesh-*")
            for subdir in sublist:
                os.chdir(subdir)
                # p1 = Popen(['lmp_mpi', '-i', 'gb4.in_final'])
                os.system("lmp_mpi -i gb4.in_final")
                os.chdir(os.pardir)
            os.chdir(os.pardir)
        return
