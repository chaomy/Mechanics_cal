#!/Users/chaomingyang/anaconda2/bin/python

import numpy as np
import os
import pickle
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS
import matplotlib.pylab as plt


class Phonon_Post(object):

    def __init__(self):
        self._qpoints = []
        self._dynamical_Matrix = []
        self._eigvc = []
        self._eigval = []
        self._band_data = None
        self._bulk = None
        self._symmetry = None
        self._force_sets = None
        self._phonon = None
        self._imaginary_q = []
        self._imaginary_band_index = []

        self.set_bulk()
        self.setup()

    def setup(self):
        phonon = Phonopy(self._bulk, supercell_matrix=[
                         [6, 0, 0], [0, 6, 0], [0, 0, 6]])
        self._phonon = phonon
        self._symmetry = self._phonon.get_symmetry()
        self._force_sets = parse_FORCE_SETS()
        phonon.get_displacement_dataset()

    def set_bulk(self):
        self._bulk = read_vasp("POSCAR")

    def append_band(self, bands, q_start, q_end):
        band = []
        for i in range(51):
            band.append(np.array(q_start) +
                        (np.array(q_end) - np.array(q_start)) / 50 * i)
        bands.append(band)

    def set_band_data(self):
        self._band_data = self.load_band_data()

    def get_band_data(self):
        return self._band_data

    def load_band_data(self):
        f = open("band.conf", 'r')
        PC = pickle.Unpickler(f)
        Band_data = PC.load()
        for i in range(len(Band_data)):
            (q, d, freq, eigv) = Band_data[i]
        return Band_data

    def generate_band_structure(self):
        bands = []
        self.append_band(bands, [0.000, 0.000, 0.500], [0.000, 0.000, 0.000])
        self.append_band(bands, [0.000, 0.000, 0.000], [0.250, -0.250, 0.250])
        self.append_band(bands, [0.250, -0.250, 0.250], [0.250, 0.250, 0.250])
        self.append_band(bands, [0.250, 0.250, 0.250], [0.000, 0.000, 0.000])
        self.append_band(bands, [0.000, 0.000, 0.000], [-0.250, 0.000, 0.500])

        #  self.append_band(bands, [-0.049, 0.049, 0.526], [0.000, 0.000, 0.000])
        #  self.append_band(bands, [0.000, 0.000, 0.000], [0.239, -0.239, 0.239])
        #  self.append_band(bands, [0.239, -0.239, 0.239], [0.190, 0.288, 0.287])
        #  self.append_band(bands, [0.190, 0.288, 0.287], [0.000, 0.000, 0.000])
        #  self.append_band(bands, [0.000, 0.000, 0.000], [-0.288, 0.049, 0.526])

        self._phonon.set_displacement_dataset(self._force_sets)
        self._phonon.produce_force_constants()

        self._phonon.set_band_structure(bands, is_eigenvectors=True,
                                        is_band_connection=True)

        q_points, distances, frequencies, eigvecs = self._phonon.get_band_structure()
        ax = self._phonon.plot_band_structure(labels=['M1', '\Gamma', 'X1', 'R', '\Gamma',
                                                      'B1', 'M2'])
        plt.savefig("band.png")
        # ax.show()

        f = open("band.conf", 'w')
        PC = pickle.Pickler(f)
        Band_data = []
        for q, d, freq, eigv in zip(q_points,
                                    distances,
                                    frequencies,
                                    eigvecs):
            for i in range(51):
                Band_data.append((q[i],
                                  d[i],
                                  freq[i],
                                  eigv[i]))
        PC.dump(Band_data)
        self._band_data = Band_data

    def find_imaginary_frequency(self):
        Band_data = self._band_data
        with open("band.info", 'w') as fid:
            for i in range(len(Band_data)):
                (q, d, freq, eigv) = Band_data[i]
                for j in range(3):
                    if freq[j] < -0.030:
                        fid.write(
                            "The index of the Imaginary Frequency is {}\n".format(j))
                        fid.write("The freq is {}".format(freq))
                        fid.write("The q point is {}".format(q))
                        fid.write("The distance is {}".format(d))
                        self._imaginary_q.append(q)
                        self._imaginary_band_index.append(j)
            fid.close()

    def generate_dos(self):
        self._phonon.set_mesh([20, 20, 20])
        self._phonon.set_thermal_properties(t_step=10,
                                            t_max=1000,
                                            t_min=0)
        self._phonon.set_total_DOS(sigma=0.1)
        # self._phonon.plot_total_DOS().show()

    def generate_modulation_of_imaginary_points(self):
        # self._phonon.set_displacement_dataset(force_sets)
        # self._phonon.produce_force_constants()

        for i in range(len(self._imaginary_q)):
            print(self._imaginary_band_index[i], self._imaginary_q[i])
            self._phonon.set_modulations([1, 1, 1],
                                         [[self._imaginary_q[i], self._imaginary_band_index[i] - 2, 1, 0]])
            self._phonon.write_modulations()
            self._phonon.write_yaml_modulations()
            os.system("cp modulation.yaml modulation-%d.yaml" % (i))
            os.system("cat modulation.yaml >>  Band_structure_report")
            os.system("cp ./MPOSCAR ./MPOSCAR-%d" % (i))

        self._phonon.write_animation(q_point=self._imaginary_q[1],
                                     anime_type='xyz',
                                     band_index=1,
                                     amplitude=10.0,
                                     num_div=100,
                                     shift=None,
                                     filename='WRe_ani')

        # for j in range(3):
        #    self._phonon.set_modulations([1,1,1],
        #                                [[self._imaginary_q[i],
        #                                  self._imaginary_band_index[i],
        #                                  3,
        #                                  0]])
        #    self._phonon.write_modulations()
        #    self._phonon.write_yaml_modulations()
        #    os.system("cat ./modulation.yaml >>  Band_structure_report")


if __name__ == '__main__':
    M = Phonon_Post()
    M.generate_band_structure()
    # M.find_imaginary_frequency()
    # M.generate_modulation_of_imaginary_points()
#    M.generate_dos()
