#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-21 16:57:45


import os
import ase.io
import numpy as np
from optparse import OptionParser
import get_data
import gn_config
import md_pot_data
import plt_drv


class cal_md_lattice(gn_config.bcc, gn_config.fcc, gn_config.hcp,
                     get_data.get_data, plt_drv.plt_drv):

    def __init__(self):
        self.pot = md_pot_data.va_pot.Nb_pbe
        get_data.get_data.__init__(self)
        plt_drv.plt_drv.__init__(self)

    def gn_temp_atoms(self):
        atoms = self.set_bcc_convention(self.pot["latbcc"], in_size=(5, 5, 5))
        self.write_lmp_config_data(atoms)
        return atoms

    def cal_lat_bcc(self):
        delt = 0.01
        npts = 100
        alat = self.pot["latbcc"]
        for i in range(-npts, npts):
            self.pot["latbcc"] = delt * i + alat
            atoms = self.set_bcc_convention(in_size=(1, 1, 1))
            self.write_lmp_config_data(atoms)
            os.system("lmp_mpi -i in.init")

    def run_lmp_lattice(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.lattice")
        os.chdir(os.pardir)

    def hcp_lattice(self):
        self.set_hcp_lattice_constant(3.1742, 5.18)
        self.set_hcp_direction()
        atoms = self.set_hcp_convention((1, 1, 1))
        ase.io.write("hcp.txt", images=atoms, format='cfg')
        self.write_lmp_config_data(atoms)

    def plt_lat(self):
        dat = np.loadtxt("out.txt")[50:150]
        da2 = np.loadtxt("ENG.log")

        idx = da2[:, 0].argsort()
        print(da2[:, 1][idx])

        self.set_111plt()
        self.ax.plot(dat[:, 0], dat[:, 3] - min(dat[:, 3]))
        self.ax.plot(da2[:, 0][idx], da2[:, 1][idx] - min(da2[:, 1]))
        self.fig.savefig("fig_engy.png", **self.figsave)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_md_lattice()

    dispatcher = {'bcc': drv.cal_lat_bcc,
                  'plt': drv.plt_lat,
                  'temp': drv.gn_temp_atoms}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

    # ./Al-Mg.eam.fs
    # a  3.18421469227388
    # c  5.18442354562984
    # Job = cal_md_lattice(in_structure='hcp',
    #                      lattice_constant=3.18,
    #                      in_element='Mg',
    #                      in_potential="./Al-Mg.eam.fs")

    #  Job.gn_temp_atoms()
    #           Nb                      W(W.eam.fs);   W.set.txt
    #  0K       3.307                   3.165200       3.1648492
    #  150K
    #  300K     3.3107                  3.165258
    #  450K
    #  600K
    #  750K
