#!/usr/local/bin/python
# encoding: utf-8
#
import   ase
import   multi_plot
import   matplotlib.pylab as plt
import   matplotlib.patches as patches
import   numpy as np
import   ase.calculators.neighborlist
import   ase.utils.geometry
import   Intro_vasp
import   os

class dd_map_dislocation(multi_plot.single_plot,
                         Intro_vasp.vasp_change_box):
    def __init__(self):
        #  self.lattice = 3.304   #  Ta
        #  self.lattice = 3.307   #  Nb
        #  self.lattice =  3.16741543    # Mo
        self.lattice = 4.05  # Al
        self.near_enighbor_cutoff = np.sqrt(3)/2. * self.lattice
        multi_plot.single_plot.__init__(self)

        self.hwidth = 2.0;
        self.hlength = 1.2;
        self.awidth = 0.04;
        Intro_vasp.vasp_change_box.__init__(self);
        return

    def read_poscar_atoms(self, in_tag):
        if in_tag == 'perf':
            atoms = ase.io.read("POSCAR_perf.vasp", format = 'vasp')
        elif in_tag == 'init':
            atoms = ase.io.read("POSCAR", format = 'vasp')
        elif in_tag == 'relaxed':
            atoms = ase.io.read("CONTCAR", format = 'vasp')
        else:
            print "###########  Error reading configuration ###########"
        return atoms

    def test_sort(self):
        poscar_atoms = ase.io.read('./lmp_init.cfg', format = 'cfg');
        perf_atoms = ase.io.read("./lmp_perf.cfg", format = 'cfg');
        atoms_tobe_changed = ase.io.read("./cfg/Al.0.cfg", format = 'cfg');

        map_list = self.map_atoms_list(atoms_tobe_changed, poscar_atoms);
        new_perf = self.sort_atoms(perf_atoms, map_list);
        ase.io.write(filename = 'new_perf.cfg',
                     images = new_perf,
                     format = 'cfg');

        return new_perf

    def draw_dd_map(self, tag=None):
        import glob
        filelist = glob.glob("./cfg/Al.110.cfg");
        if tag is None:
            tag = 'relaxed';

        ###########################################################
        ####                   read  data                      ####
        ###########################################################
        #  self.test_sort();

        atoms_perf = ase.io.read("./new_perf.cfg", format = 'cfg');

        if tag == 'relaxed':
            atoms_relax = ase.io.read(filelist[-1], format = 'cfg');
            atoms_perf.wrap(pbc=[1,1,1]);
            atoms_relax.wrap(pbc=[1,1,1]);

        elif tag == 'init':
            atoms_relax = ase.io.read(filelist[0], format = 'cfg');
            atoms_perf.wrap(pbc=[1,1,1]);
            atoms_relax.wrap(pbc=[1,1,1]);

        fig = plt.figure(figsize=(40, 15));
        ax = fig.add_subplot(111, aspect = 'equal')

        positions_perf = atoms_perf.get_positions()
        positions_relax = atoms_relax.get_positions()

        ###########################################################
        ####               original center                     ####
        ###########################################################
        #  cutoff_list = 1.385 * np.ones(len(atoms_relax));

        #  print displacementNorm
        supercell_base  = atoms_perf.get_cell()

        draw_cir = []
        plt.xlim(np.min(positions_perf[:, 0]), np.max(positions_perf[:, 0]))
        plt.ylim(np.min(positions_perf[:, 1]), np.max(positions_perf[:, 1]))

        circlesize = 0.28

        ###########################################################
        ####            two  stacking layers                  ####
        ###########################################################
        print supercell_base
        print len(atoms_perf);

        (layerList, distance) = ase.utils.geometry.get_layers(atoms_perf, [0, 0, 1], tolerance=0.01); ### the normal of z direction
        #  (layerList, distance) = ase.utils.geometry.get_layers(atoms_perf, [0, 15265, 27420], tolerance=0.1526); ### the normal of z direction

        print layerList
        print distance

        #  draw_cir.append(patches.Circle((discenter1[0], discenter1[1]),
                                        #  circlesize,
                                        #  fill= True,
                                        #  color = 'y',
                                        #  linewidth = 3.0,
                                        #  linestyle = 'solid'))

        #  draw_cir.append(patches.Circle((discenter2[0], discenter2[1]),
                                        #  circlesize,
                                        #  fill= True,
                                        #  color = 'c',
                                        #  linewidth = 3.0,
                                        #  linestyle = 'solid'))

        displacement_magnitude_list = [];
        layer1_index = [];
        layer2_index = [];

        print layerList
        for i in range(len(atoms_perf)):
            displacement_magnitude_list.append(np.linalg.norm(positions_perf[i, 2] - positions_perf[i, 2]));

            if layerList[i] == 0:
                layer1_index.append(i);
                draw_cir.append(patches.Circle((positions_perf[i, 0], positions_perf[i, 1]),
                                                circlesize,
                                                fill= False,
                                                color = 'b',
                                                linewidth = 2.5,
                                                linestyle = 'solid'))

            if layerList[i] == 1:
                layer2_index.append(i);
                draw_cir.append(patches.Circle((positions_perf[i, 0], positions_perf[i, 1]),
                                                circlesize,
                                                fill= False,
                                                color = 'g',
                                                linewidth = 2.5,
                                                linestyle = 'solid'))


        ###########################################################
        ####            darw the solute  atom                  ####
        ###########################################################
        draw_cir.append(patches.Circle((positions_perf[-1, 0], positions_perf[-1, 1]),
                                        circlesize,
                                        fill= False,
                                        color = 'm',
                                        linewidth = 4.5,
                                        linestyle = 'solid'))

        ###########################################################
        ####            draw arrows with neighbors             ####
        ###########################################################
        cutoff_list = 1.3 * np.ones(len(atoms_perf)); # nearlist neighbor

        nl = ase.calculators.neighborlist.NeighborList(cutoff_list,
                self_interaction=False, bothways=True);

        nl.build(atoms_perf);
        nl.update(atoms_perf);
        #  cutoff_dis = np.sqrt(3)/2. * 3.28;
        cutoff_dis = 1.05 * self.lattice;

        if tag == 'init':
            arrow_coeff = 0.20;
        else:
            arrow_coeff = 0.22;

        loc_scale = 1.4;

        len(layer1_index);
        halfz = 0.48 * supercell_base[2, 2];

        max_delta = 0;

        ### arrow from blue circle
        for i in layer1_index:  # draw layer 1 to others
            indices, offsets = nl.get_neighbors(i); # all nearest neighbor atoms
            print len(indices);

            for j in layer2_index:
                if (np.linalg.norm(positions_perf[i, 0:2] - positions_perf[j, 0:2]) < cutoff_dis):  # it is the neighbor

                    # z direction
                    dz1 = positions_relax[i, 2] - positions_perf[i, 2];
                    dz2 = positions_relax[j, 2] - positions_perf[j, 2];

                    delta_mag = dz2 - dz1;
                    print "delta mag", delta_mag

                    if delta_mag > halfz:
                        delta_mag -= supercell_base[2, 2];
                    if delta_mag < -halfz:
                        delta_mag += supercell_base[2, 2];

                    if delta_mag > max_delta:
                        max_delta = delta_mag;

                    # vector direction
                    delta_vector = positions_perf[j, 0:2] - positions_perf[i, 0:2];

                    lx = arrow_coeff * delta_vector[0];
                    ly = arrow_coeff * delta_vector[1];

                    dx = delta_mag * lx;
                    dy = delta_mag * ly;

                    x = positions_perf[i, 0] +  0.5 * delta_vector[0];
                    y = positions_perf[i, 1] +  0.5 * delta_vector[1];

                    if delta_mag >= 0:
                        ax.quiver(x, y, dx, dy,
                                alpha = 0.8,
                                scale= loc_scale,
                                pivot='mid',
                                color = 'k',
                                width= self.awidth,
                                units = 'inches');

        ### arrow from green circle
        for i in layer2_index:  # draw layer 1 to others
            indices, offsets = nl.get_neighbors(i); # all nearest neighbor atoms
            #  print len(indices);

            for j in layer1_index:
                if (np.linalg.norm(positions_perf[i, 0:2] - positions_perf[j, 0:2]) < cutoff_dis):  # it is the neighbor

                    # z direction
                    dz1 = positions_relax[i, 2] - positions_perf[i, 2];
                    dz2 = positions_relax[j, 2] - positions_perf[j, 2];

                    delta_mag = dz2 - dz1;

                    if delta_mag > halfz:
                        delta_mag -= supercell_base[2, 2];
                    if delta_mag < -halfz:
                        delta_mag += supercell_base[2, 2];

                    # vector direction
                    delta_vector = positions_perf[j, 0:2] - positions_perf[i, 0:2];

                    lx = arrow_coeff * delta_vector[0];
                    ly = arrow_coeff * delta_vector[1];

                    dx = delta_mag * lx;
                    dy = delta_mag * ly;

                    x = positions_perf[i, 0] +  0.5 * delta_vector[0];
                    y = positions_perf[i, 1] +  0.5 * delta_vector[1];

                    if delta_mag >= 0:
                        ax.quiver(x, y, dx, dy,
                                alpha = 0.8,
                                scale= loc_scale,
                                pivot='mid',
                                color = 'r',
                                width= self.awidth,
                                units = 'inches');

        #  ### start from red
        #  for i in layer3_index:  # draw layer 1 to others
            #  indices, offsets = nl.get_neighbors(i); # all nearest neighbor atoms
            #  #  print len(indices);

            #  for j in layer1_index:
                #  if (np.linalg.norm(positions_perf[i, 0:2] - positions_perf[j, 0:2]) < cutoff_dis):  # it is the neighbor

                    #  # z direction
                    #  dz1 = positions_relax[i, 2] - positions_perf[i, 2];
                    #  dz2 = positions_relax[j, 2] - positions_perf[j, 2];

                    #  delta_mag = dz2 - dz1;

                    #  if delta_mag >= halfz:
                        #  delta_mag -= supercell_base[2, 2];
                    #  if delta_mag <= -halfz:
                        #  delta_mag += supercell_base[2, 2];

                    #  # vector direction
                    #  delta_vector = positions_perf[j, 0:2] - positions_perf[i, 0:2];

                    #  lx = arrow_coeff * delta_vector[0];
                    #  ly = arrow_coeff * delta_vector[1];

                    #  dx = delta_mag * lx;
                    #  dy = delta_mag * ly;

                    #  x = positions_perf[i, 0] +  0.5 * delta_vector[0];
                    #  y = positions_perf[i, 1] +  0.5 * delta_vector[1];

                    #  if delta_mag >= 0:
                        #  ax.quiver(x, y, dx, dy,
                                #  alpha = 0.8,
                                #  scale= loc_scale,
                                #  pivot='mid',
                                #  color = 'k',
                                #  width= self.awidth,
                                #  units = 'inches');

            #  for k in layer2_index:
                #  if (np.linalg.norm(positions_perf[i, 0:2] - positions_perf[k, 0:2]) < cutoff_dis):  # it is the neighbor

                    #  # z direction
                    #  dz1 = positions_relax[i, 2] - positions_perf[i, 2];
                    #  dz2 = positions_relax[k, 2] - positions_perf[k, 2];

                    #  delta_mag = dz2 - dz1;

                    #  if delta_mag > halfz:
                        #  delta_mag -= supercell_base[2, 2];
                    #  if delta_mag < -halfz:
                        #  delta_mag += supercell_base[2, 2];

                    #  # vector direction
                    #  delta_vector = positions_perf[k, 0:2] - positions_perf[i, 0:2];

                    #  lx = arrow_coeff * delta_vector[0];
                    #  ly = arrow_coeff * delta_vector[1];

                    #  dx = delta_mag * lx;
                    #  dy = delta_mag * ly;

                    #  x = positions_perf[i, 0] +  0.5 * delta_vector[0];
                    #  y = positions_perf[i, 1] +  0.5 * delta_vector[1];

                    #  if delta_mag >= 0:
                        #  ax.quiver(x, y, dx, dy,
                                #  alpha = 0.8,
                                #  scale= loc_scale,
                                #  pivot='mid',
                                #  color = 'k',
                                #  width= self.awidth,
                                #  units = 'inches');

        for cir in draw_cir:
            ax.add_patch(cir)

        print len(layer1_index), len(layer2_index)
        print max_delta

        figName = "%s.png"%(tag);
        plt.savefig(figName,
                    bbox_inches='tight', pad_inches=0.01);
        #  plt.show()
        return

if __name__ == '__main__':
    Job = dd_map_dislocation()
    Job.draw_dd_map('relaxed')
    #  Job.draw_dd_map('init')
