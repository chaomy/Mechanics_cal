#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : md_pot_data.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

import numpy as np


class md_pot:
    Nb_eam = {'element': 'Nb',
              'latbcc': 3.30789893314849,
              'structure': 'bcc',
              'pair_style': 'eam/alloy',
              'pair_type': 'eam',
              'lattice': 3.30789893314849,
              'C11': 244,
              'C12': 136,
              'C44': 32,
              'file': "Nb.eam.alloy.webarchive"}

    Nb_adp = {'element': 'Nb',
              'structure': 'bcc',
              'latbcc': 3.30409,
              'pair_style': 'adp',
              'pair_type': 'adp',
              'lattice': 3.30409,
              'file': "dummy.lammps.ADP",
              'ebcc': -7.00999711903546,
              'fvac': 2.63434437728,
              'latfcc': 4.1805780840553,
              'efcc': -6.79110430035352,
              'ahcp': 2.95651917237571,
              'chcp': 4.82434648200844,
              'ehcp': -6.79091689477516,
              'C11': 258.80688,
              'C12': 128.52838,
              'C44': 20.26741,
              'B': 171.9545}

    Nb_adp_tmp = {'element': 'Nb',
                  'structure': 'bcc',
                  'latbcc': 3.31033,
                  'pair_style': 'adp',
                  'pair_type': 'adp',
                  'lattice': 3.31033,
                  'file': "dummy.lammps.ADP",
                  'ebcc': -7.00999711903546,
                  'fvac': 2.63434437728,
                  'latfcc': 4.1805780840553,
                  'efcc': -6.79110430035352,
                  'ahcp': 2.95651917237571,
                  'chcp': 4.82434648200844,
                  'ehcp': -6.79091689477516,
                  'C11': 258.80688,
                  'C12': 128.52838,
                  'C44': 20.26741,
                  'B': 171.9545}

    W_eam = {'element': 'W',
             'structure': 'bcc',
             'lattice': 3.14339,
             'latbcc': 3.14339,
             'pair_style': 'eam/alloy',
             'pair_type': 'eam',
             'file': 'w_eam4.fs'}


class dft_data:
    Nb_pbe = {'element': 'Nb',
              'lattice': 3.3224040,
              'latbcc': 3.3224040,
              'latfcc': 4.23405,
              'ahcp': 2.891,
              'chcp': 5.266,
              'C11': 246.,
              'C12': 137.,
              'C44': 20.,
              'ebcc': -10.0898432655645,
              'esingle': -3.09477080931,
              'efcc': -9.7707483945,
              'ehcp': -9.79518390985,
              'kpoints': (31, 31, 31)}

    W_pbe = {'element': 'W',
             'lattice': 3.17205,
             'latbcc': 3.17205,
             'latfcc': 4.02425,
             'ahcp': 2.7631994141,
             'chcp': 4.926472254395555,
             'ebcc': -13.0169747065,
             'ehcp': -12.50878207695,
             'efcc': -12.52477743}


class qe_pot:
    vca_W50Re50 = {'element': 'W',
                   'lattice': 3.1485,
                   'latticebohr': 5.94977,
                   'latbcc': 3.1485,
                   'mass': 183.95835,
                   'file': 'WRe.0-50.fhi.UPF'}

    vca_W75Re25 = {'element': 'W',
                   'lattice': 3.158786,
                   'latticebohr': 5.96924,
                   'mass': 184.431750,
                   'latbcc': 3.158786,
                   'file': 'WRe.0-25.fhi.UPF'}

    vca_W80Re20 = {'element': 'W',
                   'lattice':  3.16097001649,
                   'latticebohr': 5.97337,
                   'mass': 184.313400,
                   'file': 'WRe.0-20.fhi.UPF'}

    vca_W85Re15 = {'element': 'W',
                   'lattice': 3.1634677319,
                   'latticebohr': 5.97809,
                   'mass': 184.19505,
                   'file': 'WRe.0-15.fhi.UPF'}

    vca_W90Re10 = {'element': 'W',
                   'lattice': 3.16596544737,
                   'latticebohr': 5.98281,
                   'mass':  184.0767,
                   'file': 'WRe.0-10.fhi.UPF'}

    vca_W95Re05 = {'file': 'WRe.0-05.fhi.UPF'}


class unitconv:
    uengy = {'rytoeV': 13.605698066,
             'evtoJ': 1.60218e-19}
    ulength = {'AtoBohr': 1.889725988579,
               'BohrtoA': 0.529177}
    ustress = {'evA3toGpa': 160.21766208}
    ufreq = {'mevtokHz': 1. / 8.0532}


class phononpath:
    primbcc = {'\Gamma': np.array([0.0, 0.0, 0.0]),
               'H': np.array([-0.5, 0.5, 0.5]),
               'P': np.array([0.25, 0.25, 0.25]),
               'N': np.array([0.0, 0.0, 0.5])}
