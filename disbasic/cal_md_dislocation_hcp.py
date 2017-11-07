#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-13 13:25:35

import ase
from numpy import sqrt


class md_dislocation_hcp(object):

    def hcp_edge_dislocation(self):
        lata = 3.2019267694893
        latc = 5.1969105399
        atoms = ase.io.read("../MgNd_2NNmeam/Mgcfg/Mg.0.cfg")

        # cut a layers to generate free surface
        atoms = self.cut_y_normal_atoms(atoms, gp_n=2)
        atoms = self.intro_single_edge_atoms(atoms)

        # cut a layer normal the burger direction
        atoms = self.cut_x_normal_atoms(atoms,
                                        lata, 1, sqrt(3) / 4.0)
        ase.io.write("edge.cfg", atoms, "cfg")

        self.write_lmp_config_data(atoms)
        os.system("cp ./lmp_init.txt ../MgNd_2NNmeam")
        return
