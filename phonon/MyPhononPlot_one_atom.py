#!/Users/chaomingyang/anaconda2/bin/python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-25 22:13:44


from optparse import OptionParser
from itertools import cycle
from numpy import genfromtxt
import pickle as pc
import matplotlib.pylab as plt
import numpy as np
import matplotlib.ticker as ticker
import plt_drv
import md_pot_data


class MyphononPlot(plt_drv.plt_drv):

    def __init__(self):
        # self.SymbolsA = ['[0.5 0.5 0.0]','$\Gamma$',
        #        '[0.5 0.0 0.0]','[0.5 0.5 0.5]','$\Gamma$']
        # self.SymbolsA = ['N', '$\Gamma$', 'P', 'H', '$\Gamma$']
        self.SymbolsA = ['$\Gamma$', 'H', 'P', '$\Gamma$', 'N']
        # self.SymbolsB = ['[0.5 0.0 0.5]','$\Gamma$',
        #        '[0.0 0.5 0.0]','[0.0 0.5 0.5]','$\Gamma$']
        self.SymbolsB = ['N2', '$\Gamma$', 'H2', 'D1', '$\Gamma$']
        self.SymbolsC = ['X1', '$\Gamma$', 'H3', 'X2', '$\Gamma$']
        self.SymbolsD = ['M1', '$\Gamma$', 'X1', 'R1', '$\Gamma$', 'B1', 'M2']
        self.SymbolsD_one = ['N', '$\Gamma$', 'SM', 'P', '$\Gamma$', 'GP']
        self.SymbolsE = ['$\Gamma$', 'X', 'Y', '$\Sigma$', '$\Gamma$',
                         'Z', '$\Sigma_1$', 'N', 'P', '$ Y_1$', 'Z', 'X', 'P']
        self.SymbolsF = ['$\Gamma$', 'X', 'L', 'T', 'W', 'R', '$ X_1$', 'Z',
                         '$\Gamma$', 'Y', 'S', 'W', '$ L_1$', 'Y', '$ Y_1$', 'Z']

        if options.mpath == 'm1':
            self.savename = 'Ph_M1.png'
            self.used_symbol = self.SymbolsD_one

        elif options.mpath == "TP":
            self.savename = 'Ph_TP.png'
            self.used_symbol = self.SymbolsE

        elif options.mpath == "OP":
            self.savename = 'Ph_OP.png'
            self.used_symbol = self.SymbolsF

        elif options.mpath == 'bcc':
            self.savename = 'bcc.png'
            self.used_symbol = self.SymbolsA

        self.ToThz = 1. / 33.3
        # self.mfigsize = (10, 6.5)

        plt_drv.plt_drv.__init__(self)
        self.set_keys()

    def set_pltkargs(self):
        keyslist = [{'linestyle': '-',
                     'color': next(self.coloriter), 'linewidth': 4},
                    {'linestyle': '--',
                     'color': next(self.coloriter), 'linewidth': 4},
                    {'linestyle': ':',
                     'color': next(self.coloriter), 'linewidth': 6}]
        return keyslist

    def plot_band_setting(self):
        self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
        self.ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))

    ##############################################################
    # plot the frequency (interface to phonopy)
    ##############################################################
    def read_config(self, config):
        banddata = np.load("{}.npy".format(config))

        freqs = []
        distances = []
        primbcc = md_pot_data.phononpath.primbcc
        sympnts = []

        for i in range(len(banddata)):
            (q, d, freq, eigv) = banddata[i]
            if (q == primbcc['\Gamma']).all():
                sympnts.append(d)

            elif (q == primbcc['H']).all():
                sympnts.append(d)

            elif (q == primbcc['P']).all():
                sympnts.append(d)

            elif (q == primbcc['N']).all():
                sympnts.append(d)

            distances.append(d)
            freqs.append(freq)

        distances = np.array(distances)
        freqs = np.array(freqs)
        return (distances, freqs, sympnts)

    def plot_freq(self, **kwargs):
        if 'config' in kwargs.keys():
            (distances, freqs, sympnts) = self.read_config(kwargs['config'])
        else:
            (distances, freqs, sympnts) = \
                kwargs['distances'], kwargs['freqs'], kwargs['sympnts']

        if 'ax' not in kwargs.keys():
            self.set_111plt()
            ax = self.ax
            ax.set_yticks(np.arange(0, 7.0, 1.5))
            # ax.set_xlim([0.0, 1.0])
            ax = self.add_vline(self.ax, sympnts)
        else:
            ax = kwargs.get('ax')

        if 'pltkeys' in kwargs.keys():
            pltkeys = kwargs['pltkeys']
        else:
            self.set_phonon_keys()
            pltkeys = self.pltkwargs
        plt.ylim(0, 7.0)
        ax.plot(distances, freqs[:, 0], label=kwargs['label'], **pltkeys)
        # ax.plot(distances, freqs[:, 0], **pltkeys)
        ax.plot(distances, freqs[:, 1], **pltkeys)
        ax.plot(distances, freqs[:, 2], **pltkeys)
        return ax

    def plot_exp(self, **kwargs):
        nodes = np.array([0.0, 0.3010743115163611, 0.5618123137164394,
                          0.8225503159165176, 1.0354420032308085])
        GH01 = genfromtxt('GH_01.csv', delimiter=',')
        GH02 = genfromtxt('GH_02.csv', delimiter=',')
        print(GH01.shape)
        GH01[:, 0] *= nodes[1]
        GH02[:, 0] *= nodes[1]
        kwargs['ax'].plot(GH01[:, 0], GH01[:, 1],
                          label='Experiment', **kwargs['pltkeys'])
        kwargs['ax'].plot(GH02[:, 0], GH02[:, 1], **kwargs['pltkeys'])

        HPG01 = genfromtxt('HPG_01.csv', delimiter=',')
        HPG02 = genfromtxt('HPG_02.csv', delimiter=',')

        HPG01[:, 0] = HPG01[:, 0] * (nodes[3] - nodes[1]) + nodes[1]
        HPG02[:, 0] = HPG02[:, 0] * (nodes[3] - nodes[1]) + nodes[1]

        kwargs['ax'].plot(HPG01[:, 0], HPG01[:, 1], **kwargs['pltkeys'])
        kwargs['ax'].plot(HPG02[:, 0], HPG02[:, 1], **kwargs['pltkeys'])

        GN01 = genfromtxt('GN_01.csv', delimiter=',')
        GN02 = genfromtxt('GN_02.csv', delimiter=',')
        GN03 = genfromtxt('GN_03.csv', delimiter=',')

        GN01[:, 0] = GN01[:, 0] * (nodes[4] - nodes[3]) + nodes[3]
        GN02[:, 0] = GN02[:, 0] * (nodes[4] - nodes[3]) + nodes[3]
        GN03[:, 0] = GN03[:, 0] * (nodes[4] - nodes[3]) + nodes[3]

        kwargs['ax'].plot(GN01[:, 0], GN01[:, 1], **kwargs['pltkeys'])
        kwargs['ax'].plot(GN02[:, 0], GN02[:, 1], **kwargs['pltkeys'])
        kwargs['ax'].plot(GN03[:, 0], GN03[:, 1], **kwargs['pltkeys'])
        return kwargs['ax']

    def plot_single(self, config):
        inargs = {'config': config}
        ax = self.plot_freq(**inargs)
        plt.yticks(size=self.myfontsize)
        plt.xticks(size=self.myfontsize)
        plt.savefig("mphonon.png", **self.figsave)

    def plot_freq_cmp(self):
        configlist = ['band.conf.vasp', 'band.conf.lmp']
        self.set_111plt()
        self.ax.set_yticks(np.arange(0, 7.0, 1.5))
        self.ax.set_xlim([0.0, 1.0])
        keyslist = self.set_pltkargs()
        self.set_phonon_keys()
        for config, lab, cnt in zip(configlist,
                                    ['PAW-PBE', 'MEAMS'],
                                    range(len(configlist))):
            (distances, freqs, sympnts) = self.read_config(config)
            inargs = {'distances': distances,
                      'freqs': freqs,
                      'sympnts': sympnts,
                      'ax': self.ax,
                      'pltkeys': keyslist[cnt],
                      'label': lab}
            self.ax = self.plot_freq(**inargs)
        self.ax = self.add_vline(self.ax, sympnts)
        self.ax = self.add_phonon_labels(self.ax, sympnts)

        # plot experiment
        self.ax = self.plot_exp(**{'ax': self.ax, 'pltkeys': keyslist[2]})
        self.add_y_labels(cycle(['Frequency [THz]']), self.ax)
        self.add_legends(self.ax)
        self.set_tick_size(self.ax)
        plt.savefig("mphonon.png", **self.figsave)

    def add_phonon_labels(self, ax, sympnts):
        print(set(sympnts).intersection())
        pnts = []
        for point in set(sympnts).intersection():
            pnts.append(point)
        pnts = np.sort(np.array(pnts))
        plt.xticks(pnts, self.SymbolsA, size=self.myfontsize)
        return ax

    def add_vline(self, ax, sympnts):
        for point in set(sympnts).intersection():
            ax.axvline(x=point, alpha=0.5,
                       linestyle='-.', linewidth=2.0, color='k')
        return ax

    ##############################################################
    # plot the smallest frequency  (interface to phonopy)
    ##############################################################
    def plot_smallest_frequency(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))

        for i in range(len(options.infile)):
            print(options.infile[i])
            fid = open(options.infile[i])
            mylabel = "$\epsilon_{11}$ = %s" % (options.infile[i][-5:-1])
            A = pc.Unpickler(fid)

            [symbols, sympnts, distances, Freqs] = A.load()
            print(distances)
            print(len(distances))
            print(sympnts)

            if i == 0:
                for point in sympnts:
                    plt.axvline(x=point,
                                linestyle=':', linewidth=2.0, color='k')

            for j in range(len(distances)):
                Freqs[j] = Freqs[j]
                if j == 1:
                    plt.plot(distances[j],
                             Freqs[j],
                             color=self.color[i], linestyle='-',
                             linewidth=2, marker='o',
                             markersize=3, label=mylabel)
                else:
                    plt.plot(distances[j],
                             Freqs[j],
                             color=self.color[i], linestyle='-',
                             linewidth=2, marker='o', markersize=3)

            plt.xlim(0, distances[-1][-1])
            plt.ylim(-1.0, 4.0)

            if symbols:
                plt.xticks(sympnts, self.used_symbol,
                           size=self.myfontsize)
            else:
                plt.xticks(sympnts, [''] * len(sympnts),
                           size=self.myfontsize)

        plt.yticks(size=self.myfontsize)
        plt.xticks(size=self.myfontsize)

        plt.legend(fontsize=self.myfontsize - 2, mode="expand",
                   borderaxespad=0)

        # plt.legend(bbox_to_anchor=(0., 1.00, 1., .100), loc=3,
        #            ncol=3, mode="expand", borderaxespad=0,
        #            fontsize=self.myfontsize - 2)

        plt.axhline(y=0, linestyle='--', linewidth=2.0,
                    color='k')
        plt.ylabel('Frequency [THz]', {'fontsize': self.myfontsize})
        plt.xlabel('Wave vector', {'fontsize': self.myfontsize})
        plt.savefig(self.savename,
                    bbox_inches='tight', pad_inches=0.02)
        plt.show()


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--infile", action="append")
    parser.add_option("-s", "--outfile", action="append")
    parser.add_option("-p", "--path", type="string",
                      dest="mpath", default='bcc')
    parser.add_option("-t", "--type", action="store",
                      type="string", dest="mtype")

    (options, args) = parser.parse_args()
    drv = MyphononPlot()
    if options.mtype.lower() == 'vaspsmall':
        drv.plot_smallest_frequency()

    elif options.mtype.lower() == 'vasp':
        drv.plot_freq('band.conf.vasp')

    elif options.mtype.lower() == 'lmp':
        drv.plot_single('band.conf.lmp')

    elif options.mtype.lower() == 'cmp':
        drv.plot_freq_cmp()

    elif options.mtype.lower() == 'exp':
        drv.plot_exp()

    ##############################################################
    # check ~/bin/Phonon/MyPhononPlot_old_version.py
    ##############################################################
    elif options.mtype == 'qe':
        drv.plot_qe_freq()

    elif options.mtype == 'qem' or options.mtype == 'mqe':
        drv.plot_qe_many()
