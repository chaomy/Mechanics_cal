#!/usr/local/bin/python
# encoding: utf-8
#
import numpy as np
import os, shutil
try:
    import  get_data
    import  gn_pbs
    import  gn_dd_ctrl
    import  dis_mob_function

except ImportError:
    print "error during import"

class dd_tensile(gn_dd_ctrl.gn_dd_ctrl,
                 gn_pbs.gn_pbs,
                 get_data.get_data,
                 dis_mob_function.fit_mob_function):

    def __init__(self):
        self.lattice_constant = 3.16;

        gn_pbs.gn_pbs.__init__(self);
        get_data.get_data.__init__(self);
        dis_mob_function.fit_mob_function.__init__(self, self.lattice_constant);

        root_dir = os.getcwd();
        return;

    def dd_mob_prepare(self):
        (tau, templist, edge_moblist, screw_moblist) = self.cal_dis_mob();

        print templist;
        print edge_moblist;
        print screw_moblist;
        fm_table = "W_miu160_niu2p79";
        config = "glide.data";

        ### set pbs ####
        self.set_nnodes(1);
        self.set_wall_time(140);
        self.set_ppn(8);
        self.set_main_job("mpirun -n 8 paradis glide.ctrl")

        for temp, edge_mob, screw_mob in zip(templist, edge_moblist,
                screw_moblist):
            file_dir = "dir-%d"%(temp);
            if not os.path.isdir(file_dir):
                os.mkdir(file_dir);

            shutil.copy(fm_table, file_dir);
            shutil.copy(config, file_dir);

            os.chdir(file_dir);
            ### write pbs ####
            self.set_job_title("tau1400-1e4-%s"%(file_dir));
            self.write_pbs()

            self.write_ctrl(screw_mob, edge_mob);

            os.chdir(self.root_dir);
        return

    def loop_sub(self):
        import glob
        dir_list = glob.glob("dir-*");

        for the_dir in dir_list:
            os.chdir(the_dir);

            os.system("qsub va.pbs");

            os.chdir(self.root_dir);
        return

if __name__ == '__main__':
    job = dd_tensile();
    job.dd_mob_prepare();
    #  job.loop_sub();
