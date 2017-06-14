#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : test.py
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

try:
    import numpy as np

except ImportError:
    print("error during import")


def myfunct(**kwargs):
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            if key is 'cat': 
                print "%s == %s" % (key, value)




        print kwargs.keys()
        print kwargs.values()


import gn_pbs
def test_prop():
    drv = gn_pbs.gn_pbs()
    print drv.nnodes
    return

if __name__ == '__main__':
    myfunct(**{'cat':1, 'dog':2})
    test_prop()
