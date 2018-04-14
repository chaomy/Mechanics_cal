#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-13 20:48:15


from optparse import OptionParser
from gn_config import gnStructure
from copy import deepcopy
from plt_drv import plt_drv
from os import system
from . import cal_va_gsf
import ase.lattice.cubic as cubic
import ase.lattice.orthorhombic as otho
import numpy as np
import gn_pbs
import get_data
import ase.io
import os
import md_pot_data
from numpy import sqrt


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 1. / 3., 1. / 2.],
                     [1 / 2., 5. / 6., 1. / 2.]]
othoHCP = othoHCPFractory()

# 1 eV = 1.60218e-19 J
coeff = 16.0218
area = 17.836829147070748


class cal_hcp_gsf(get_data.get_data,
                  gn_pbs.gn_pbs,
                  gnStructure,
                  plt_drv,
                  cal_va_gsf.cal_va_gsf):

    def __init__(self):
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gnStructure.__init__(self)
        plt_drv.__init__(self)
        cal_va_gsf.cal_va_gsf.__init__(self)

        # self.pot = md_pot_data.md_pot.mg_kim
        self.pot = md_pot_data.va_pot.Mg_pbe
        self.disps = np.arange(0.0, 0.5, 0.05)

    def gn_displacement(self, atoms,
                        displacement_vector):
        positions = atoms.get_positions()
        atom_num = len(positions)
        displacement = deepcopy(positions)
        cut = 0.5 * np.max(positions[:, 2])
        for i in range(atom_num):
            if positions[i, 2] < cut:
                displacement[i] = [0, 0, 0]
            else:
                displacement[i] = displacement_vector
        return displacement

    def set_pbs(self, mdir, opt='qe'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("{}".format(mdir))
        self.set_wall_time(100)
        self.set_main_job("""mpirun vasp > vasp.log""")
        self.write_pbs(od=True)

    def buildhcp(self):
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'],
                                         self.pot['ahcp'] * sqrt(3.),
                                         self.pot['chcp']),
                        size=(1, 1, 18),
                        symbol=self.pot['element'])
        self.write_lmp_config_data(atoms)
        return atoms

    def clc_data(self, ptype="lmp"):
        disps = self.disps
        npts = len(disps)
        raw = np.ndarray([npts, 3])
        area = self.pot["ahcp"] * self.pot["ahcp"] * sqrt(3.)
        for i, disp in zip(list(range(npts)), disps):
            mdir = 'dir_{}_{:4.3f}'.format(i, disp)
            os.chdir(mdir)

            raw[i, 0] = disp
            if ptype == "lmp":
                raw[i, 1] = np.loadtxt("out.txt")
            elif ptype == "vasp":
                (energy, vol, atoms) = self.vasp_energy_stress_vol_quick()
                raw[i, 1] = energy
                raw[i, 2] = energy / area
            os.chdir(os.pardir)
        print(raw)
        np.savetxt("save.txt", raw)

    def trans(self):
        disps = self.disps
        npts = len(disps)
        for i, disp in zip(list(range(npts)), disps):
            mdir = 'dir_{}_{:4.3f}'.format(i, disp)
            PTH = "$FLUX:/scratch/qiliang_flux/chaomy/VA/Mg/gsfMgVa1-120"
            system("scp {}/{}/OUTCAR {}".format(PTH, mdir, mdir))
            system("scp {}/{}/CONTCAR {}".format(PTH, mdir, mdir))

    def clc_plt(self):
        raw = np.loadtxt("save.txt")
        area = self.pot["ahcp"] * self.pot["ahcp"] * sqrt(3.)
        raw[:, 1] = (raw[:, 1] - raw[0, 1]) / area
        self.set_111plt()
        self.set_keys()
        self.ax.plot(raw[:, 0], raw[:, 1], **next(self.keysiter))
        self.fig.savefig("fig_gsf.png", **self.figsave)
        inev = max(raw[:, 1])
        # Nb->0.05 eV / A ^ 2
        print("ugsf = ",  inev, "ev / a^2;", inev * coeff, "J/m^2")

    def cal_stacking(self, opt="pre"):
        atoms = self.buildhcp()
        perf_cells = deepcopy(atoms.get_cell())
        self.write_lmp_config_data(atoms)
        disps = self.disps
        npts = len(disps)
        for i, disp in zip(list(range(npts)), disps):
            mdir = 'dir_{}_{:4.3f}'.format(i, disp)
            self.mymkdir(mdir)
            os.chdir(mdir)
            if opt in ["run"]:
                os.system("mpirun -n 4 lmp_mpi -i in.init")
            else:
                disp_vector = [disp, disp, 0.0]  # 1 - 120
                # disp_vector = [0.0, disp, 0.0]  # 10-10

                disp_matrix_direct = self.gn_displacement(
                    atoms.copy(), disp_vector)
                disp_matrix = deepcopy(disp_matrix_direct)
                disp_matrix[:, 0] = disp_matrix_direct[
                    :, 0] * perf_cells[0, 0]
                disp_matrix[:, 1] = disp_matrix_direct[
                    :, 1] * perf_cells[1, 1]

                local_atoms = atoms.copy()
                local_atoms.translate(disp_matrix)

                self.set_pbs(mdir)
                # self.write_lmp_config_data(local_atoms)

                ase.io.write("POSCAR", images=local_atoms, format='vasp')
                self.add_selective_dynamics(len(disp_matrix))
                os.system("cp POSCAR ../lmp.{:03d}".format(i))
                os.system("cp ../INCAR ../KPOINTS  ../POTCAR .")

            os.chdir(os.pardir)
        print(("area = ", self.pot["ahcp"] * self.pot["ahcp"] * sqrt(3.)))


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_hcp_gsf()
    dispatcher = {'prep': drv.buildhcp,
                  'gsf': drv.cal_stacking,
                  'clc': drv.clc_data,
                  'plt': drv.clc_plt,
                  'trans': drv.trans}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
