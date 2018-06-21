#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-06 15:42:15
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-03 16:00:39


from optparse import OptionParser
import numpy as np
import ase
import gn_config
import get_data
import glob
import os
import md_pot_data
import gn_pbs
import plt_drv


class cal_md_twin(gn_config.gnStructure, get_data.get_data,
                  plt_drv.plt_drv, gn_pbs.gn_pbs):

    def __init__(self):
        get_data.get_data.__init__(self)
        self.pot = self.load_data("../BASICS/pot.dat")
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)

    def perfbcc(self, layer=8):
        atoms = self.set_bcc_convention([[1, 1, 1], [0, -1, 1],
                                         [2, -1, -1]], (1, 1, layer))
        return atoms

    def set_pbs(self, mdir):
        self.set_pbs_type('va')
        self.set_wall_time(40)
        self.set_job_title(mdir)
        self.set_nnodes(4)
        self.set_ppn(12)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=False)

    def prepare_vasp_inputs(self, mdir):
        self.set_pbs(mdir)
        os.system("mv va.pbs {}".format(mdir))
        os.system("cp KPOINTS {}".format(mdir))
        os.system("cp INCAR {}".format(mdir))
        os.system("cp POTCAR {}".format(mdir))

    def bcc211(self):
        # layer = 3 
        layer = 4
        total = 12
        atoms = self.perfbcc(total)
        la = self.pot["latbcc"]
        uz = la * np.sqrt(6) / 6.
        ux = la * np.sqrt(3) / 6.
        pos = atoms.get_positions()
        disp = np.zeros(pos.shape)

        # add vacumm
        cell = atoms.get_cell()
        # print("vacumm", uz * 10)
        # cell[2, 2] += uz * 12
        # atoms.set_cell(cell)

        opt = 'md'
        startshif = 4 

        for shif in range(layer * 6):
            atoms_shift = atoms.copy()
            for i in range((total - startshif) * 6):
                if (i < shif):
                    disp[startshif * 6 + i, 0] = i * ux
                else:
                    disp[startshif * 6 + i, 0] = shif * ux
            atoms_shift.translate(disp)
            atoms_shift.set_positions(atoms_shift.get_positions())
            # atoms_shift.set_positions(atoms_shift.get_positions() +
            #                           np.array([0., 0., 6 * uz]))
            if opt in ['md']:
                self.write_lmp_config_data(atoms_shift,
                                           "pos_{:03}".format(shif))
            elif opt in ['va']:
                mdir = 'dir_{:03d}'.format(shif)
                self.mymkdir(mdir)
                ase.io.write("{}/POSCAR".format(mdir), atoms_shift, 'vasp')
                self.prepare_vasp_inputs(mdir)

    def run(self):
        fls = glob.glob("pos_*")
        print(fls)
        os.system("cp {} lmp_init.txt".format(fls[0]))
        os.system("lmp_mpi -i in.init > logi")
        os.system("mv out.txt outi.txt")
        os.system("cp {} lmp_init.txt".format(fls[-1]))
        os.system("lmp_mpi -i in.init > logf")
        os.system("mv out.txt outf.txt")

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

    def twin_vasp(self):
        coeff = 16.021766208
        area = 2.8771910028406178 * 4.6984332329907597
        twin = -419.44619583
        notwin = -419.90847454
        twindft = 0.5 * (twin - notwin) / area
        # 0.27394446746637624
        print("DFT", 1e3 * twindft, "meV / A^2",
              1e3 * twindft * coeff, "mJ / m^2")
        atoms = self.perfbcc(15)
        cell = atoms.get_cell()
        area = cell[0, 0] * cell[1, 1]
        twin = np.loadtxt("outf.txt")
        notwin = np.loadtxt("outi.txt")
        twine = 0.5 * (twin - notwin) / area
        print("MEAMS", 1e3 * twine, "meV / A^2",
              1e3 * twine * coeff, "mJ / m^2")
        print("ERROR", (twine - twindft) / twindft)

    # twin
    # DFT   twin energy : -419.44619583
    # DFT   no tiwn : -419.90847454

    def auto(self):
        self.bcc211()
        self.run()
        self.twin_vasp()

    def transpose_dump_to_poscar(self):
        atoms = ase.io.read("dump.rst", format='lammps-dump') 
        ase.io.write("poscar.rst", images=atoms, format='vasp')

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
                  'plt': drv.plt,
                  'vasp': drv.twin_vasp,
                  'auto': drv.auto,
                  'trans': drv.transpose_dump_to_poscar}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
