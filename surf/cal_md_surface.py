#!/ usr / bin / env python
#- * - coding : utf - 8 - * -
#@Author : chaomy
#@Date : 2017 - 06 - 25 14 : 28 : 58
#@Last Modified by : chaomy
#@Last Modified time : 2018 - 03 - 20 13 : 34 : 11

from multiprocessing import Pool
from optparse import OptionParser
import os
import gn_config
import get_data
import gsf_data
import gn_lmp_infile
import md_pot_data


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_surface.run_lmp_minimize(*arg, **kwarg)

mp = {'100': 'x100z100',
      '110': 'x110z110',
      '111': 'x112z111'}


class cal_md_surface(gn_config.gnStructure,
                     get_data.get_data,
                     gn_lmp_infile.gn_md_infile):

    def __init__(self, pot=None):
        self.pot = pot
        get_data.get_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)

    def gn_surface_atoms(self, tg):
        atoms = self.set_bcc_convention(
            gsf_data.gsfbase[mp[tg]], gsf_data.gsfsize[mp[tg]])
        for i in range(gsf_data.gsfpopn[mp[tg]]):
            atoms.pop()
        return atoms

    def gn_bulk_atoms(self, tg):
        atoms = self.set_bcc_convention(
            gsf_data.gsfbase[mp[tg]], gsf_data.bulksize[mp[tg]])
        return atoms

    def run_lmp_minimize(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.minimize")
        os.chdir(os.pardir)

    def loop_run_surface(self):
        dir_list = []
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            dir_list.append('dir-surf-%s' % (loop_list[i]))
            dir_list.append('dir-bulk-%s' % (loop_list[i]))
        for mdir in dir_list:
            self.run_lmp_minimize(mdir)

    def multi_thread_surface(self):
        dir_list = []
        loop_list = ['100', '110', '111']
        for i in range(len(loop_list)):
            dir_list.append('dir-surf-%s' % (loop_list[i]))
            dir_list.append('dir-bulk-%s' % (loop_list[i]))
        num_threads = len(dir_list)
        pool = Pool(processes=num_threads)
        pool.map(unwrap_self_run_lammps,
                 list(zip([self] * num_threads, dir_list)))

    def loop_cal_surface_data(self):
        surface_energy_list = []
        for ee in ['100', '110', '111']:
            surface_energy_list.append(self.cal_surface_energy(ee))
        return surface_energy_list

    def loop_prepare_surface(self):
        for ee in ['100', '110', '111']:
            self.prepare_md_surface(ee)
            self.prepare_md_bulk(ee)

    def prepare_md_surface(self, tg):
        mdir = 'dir-surf-{}'.format(tg)
        self.mymkdir(mdir)
        self.write_lmp_config_data(self.gn_surface_atoms(tg))
        os.system("cp lmp_init.txt in.minimize {}".format(mdir))

    def prepare_md_bulk(self, tg):
        mdir = 'dir-bulk-{}'.format(tg)
        self.mymkdir(mdir)
        self.write_lmp_config_data(self.gn_bulk_atoms(tg))
        os.system("cp lmp_init.txt in.minimize {}".format(mdir))

    def cal_surface_energy(self, tg):
        os.chdir('dir-bulk-{}'.format(tg))
        energy_b = self.md_get_final_energy()
        super_cell = self.md_get_cell()
        xy_area = self.cal_poscar_xy_area(super_cell)
        os.chdir(os.pardir)

        os.chdir('dir-surf-{}'.format(tg))
        energy_s = self.md_get_final_energy()
        os.chdir(os.pardir)

        surface_e = 0.5 * (energy_s - energy_b) / xy_area
        with open("surface_energy.dat", 'a') as fid:
            fid.write("surface energy %s %6.4f ev/A \n" % (tg, surface_e))
            fid.write("surface energy %s %6.4f J/m \n" % (tg,
                                                          surface_e * self.ev_angstrom2_to_j_m2))
        return (surface_e * self.ev_angstrom2_to_j_m2)

    def wrap_run(self):
        self.loop_prepare_surface()
        self.multi_thread_surface()
        print(self.loop_cal_surface_data())


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_md_surface()
    dispatcher = {'run': drv.wrap_run}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
