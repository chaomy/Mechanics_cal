#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-06 15:42:15
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-13 20:00:26


from optparse import OptionParser
import numpy as np
import ase
import gn_config
import get_data
import glob
import os
import md_pot_data
import plt_drv


class cal_md_twin(gn_config.gnStructure,
                  get_data.get_data,
                  plt_drv.plt_drv):

    def __init__(self):
        get_data.get_data.__init__(self)
        self.pot = self.load_data("../BASICS_MO/pot.dat")
        plt_drv.plt_drv.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)

    def perfbcc(self, layer=8):
        # six atoms per cell
        atoms = self.set_bcc_convention([[1, 1, 1], [0, -1, 1],
                                         [2, -1, -1]], (1, 1, layer))
        return atoms

    def bcc211(self):
        layer = 4
        atoms = self.perfbcc(2 * layer)
        la = self.pot["latbcc"]
        uz = la * np.sqrt(6) / 6.
        ux = la * np.sqrt(3) / 6.
        pos = atoms.get_positions()
        disp = np.zeros(pos.shape)

        # lambda = INF
        # for i in range(layer * 6):
        #     print(6 * layer + i * uz, i * ux)
        #     disp[layer * 6 + i, 0] = i * ux

        cell = atoms.get_cell()
        cell[2, 2] += layer * 6
        atoms.set_cell(cell)

        for shif in range(layer * 6 - 6):
            atoms_shift = atoms.copy()
            for i in range(layer * 6):
                if (i < shif):
                    disp[layer * 6 + i, 0] = i * ux
                else:
                    disp[layer * 6 + i, 0] = shif * ux
            atoms_shift.translate(disp)
            self.write_lmp_config_data(atoms_shift, "pos_{:03}".format(shif))

    def run(self):
        fls = glob.glob("pos_*")
        data = np.zeros(len(fls))
        for i in range(len(fls)):
            os.system("cp {} lmp_init.txt".format(fls[i]))
            os.system("lmp_mpi -i in.init")
            data[i] = np.loadtxt('out.txt')
        np.savetxt("twin.dat", data, fmt="%.12f")

    def plt(self):
        atoms = self.perfbcc(2)
        cell = atoms.get_cell()
        area = cell[0, 0] * cell[1, 1]
        self.set_111plt()
        self.set_keys()
        data = np.loadtxt("twin.dat")
        data = data - np.min(data)
        data *= 16.021766208
        data *= 1e3
        # data /= 2.0
        data /= area
        self.ax.plot(data, **next(self.keysiter))
        self.fig.savefig("fig_twin.png", **self.figsave)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_md_twin()
    dispatcher = {'bcc211': drv.bcc211,
                  'run': drv.run,
                  'plt': drv.plt}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
