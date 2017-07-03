#!/usr/bin/env python
# encoding: utf-8

import os
import shutil
import glob
import md_pot_data
import gn_config
import gn_pbs
import gn_qe_inputs
import get_data
from optparse import OptionParser


class cal_qe_phonon(gn_config.bcc,
                    gn_pbs.gn_pbs,
                    gn_qe_inputs.gn_qe_infile,
                    get_data.get_data):

    def __init__(self):
        self.pot = md_pot_data.qe_pot.vca_W75Re25
        gn_pbs.gn_pbs.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        get_data.get_data.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        self._scf_exe = 'mpirun pw.x < scf.in > scf.out'
        self._ph_exe = 'mpirun ph.x < ph.in > ph.out'
        self._q2r_exe = 'q2r.x < q2r.in > q2r.out'
        self._matdyn_exe = 'matdyn.x < matdyn_disp.in > matdyn_disp.out'
        self._nq1, self._nq2, self._nq3 = 6, 6, 6
        return

    def set_pbs(self, dirname='qe', opt='ph'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("{}_{}".format(opt, dirname))
        self.set_wall_time(30)
        self.set_main_job("""
        mpirun pw.x <  qe.acc.in  > qe.acc.out
                          """.format(opt))
        self.write_pbs(od=True)
        return

    def write_phinfile(self, opt='scratch'):
        with open('ph.in', 'w') as fid:
            fid.write("""
&inputph
 outdir = './results/',
 prefix = 'qe',
 ldisp  = .true.
 fildyn ='dynmat'
 nq1 =%d, nq2=%d, nq3=%d,
 tr2_ph=1.0d-13,
 /
            """ % (self._nq1, self._nq2, self._nq3))
            fid.close()
        return

    def write_q2r_infile(self):
        with open("q2r.in", 'w') as fid:
            fid.write("""
&input
  fildyn = 'dynmat'
  zasr   = 'simple'
  flfrc  = '%d%d%d.fc'
/
        """ % (self._nq1,
               self._nq2,
               self._nq3))
        return

    def write_matdyn(self):
        with open('matdyn_disp.in', 'w') as fid:
            fid.write("""
  &input
  asr='simple'
  amass(1) = {}
  flfrc = '{}{}{}.fc'
  flfrq = 'qe.freq'
  q_in_band_form = .True.
/
10
0.50000  0.5000    0.00000   10
0.00000  0.00000   0.00000   10
0.50000  0.00000   0.00000   10
0.5      0.5       0.5       10
0.00000  0.00000   0.00000   10
0.00000  0.50000   0.00000   10
0.5      0.5       0.5       10
0.00000  0.50000   0.50000   10
0.00000  0.00000   0.00000   10
0.50000  0.00000   0.50000   10
        """ % (self.pot['mass'],
               self._nq1, self._nq2, self._nq3))
        return

    def setup_scf_acc(self, opt='one'):
        if opt == 'reloop':
            dirlist = glob.glob('dir-*')
            for mdir in dirlist:
                os.chdir(mdir)
                self.setup_scf_acc('one')
                self.set_pbs(mdir, 'scf')
                os.chdir(os.pardir)
        else:
            self.gn_qe_restart()
        return

    def loop_phonon(self):
        rootdir = os.getcwd()
        infilelist = glob.glob("new-WRe*")
        for infile in infilelist:
            dirname = "Dir-%s" % (infile)
            os.mkdir(dirname)
            shutil.copy(self._Potentialist[self._Num], dirname)

            os.system("cp %s %s/WRe_scf.in" % (infile, dirname))
            os.chdir(dirname)

            self.write_phinfile()
            self.write_q2r_infile()
            self.write_matdyn()

            os.system("%s" % (self._scf_exe))
            os.system("%s" % (self._ph_exe))
            os.system("%s" % (self._q2r_exe))
            os.system("%s" % (self._matdyn_exe))

            self.clean_qe_results()
            os.chdir(rootdir)
        return

    def clean_qe_results(self):
        if os.path.isfile("WRe.freq"):
            shutil.rmtree("results")
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype",
                      default="curv")
    (options, args) = parser.parse_args()
    drv = cal_qe_phonon()
    if options.mtype.lower() in ['scf_re', 'scf_reloop']:
        opt = options.mtype.lower().split('_')[-1]
        drv.setup_scf_acc(opt)
