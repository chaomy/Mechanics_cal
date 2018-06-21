#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-11-01 01:54:26
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-09 11:46:31

from numpy import sqrt


class md_gb_lmp(object):

    def __init__(self):
        self.tilt = {'0001': self.hcp_til_mtx_z0001,
                     '1100': self.hcp_til_mtx_z1100,
                     '1120': self.hcp_til_mtx_z1120}

        self.zunit = {'1100': sqrt(3) * self.pot['ahcp'],
                      '1120': self.pot['ahcp'],
                      '0001': self.pot['chcp']}

    def build_hcp_gb_lmp(self, ang=15, xdim=10, thk=(40, 60), disp=(0, 0, 0), key='1100'):
        fid = open("in.gb", 'w')
        pmat = self.tilt[key](ang)
        nmat = self.tilt[key](-ang)
        fid.write(inmdhcp % (self.pot['ahcp'], self.pot['chcp'],
                             xdim, self.zunit[key],
                             thk[1], thk[1], thk[0], thk[0],
                             pmat[0, 0], pmat[0, 1], pmat[0, 2],
                             pmat[1, 0], pmat[1, 1], pmat[1, 2],
                             pmat[2, 0], pmat[2, 1], pmat[2, 2],
                             self.bs,
                             nmat[0, 0], nmat[0, 1], nmat[0, 2],
                             nmat[1, 0], nmat[1, 1], nmat[1, 2],
                             nmat[2, 0], nmat[2, 1], nmat[2, 2],
                             self.bs,
                             self.pot['pair_style'],
                             self.pot['file'], disp[0], disp[2], self.pot['ehcp']))
        fid.close()

inhcp = """
units metal
atom_style atomic
boundary  p  p  p 

variable  alat  equal    %f
variable  clat  equal    %f
variable  ysize equal    15 
variable  zsize equal    4 
variable  xdim  equal    %f
variable  zdim  equal    %f*${zsize}
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

variable ydimmin equal    	 -%d
variable ydimmax equal   	  %d
variable yfreedimmin equal   -%d
variable yfreedimmax equal    %d
variable ytop    equal        ${ydimmax}+15
variable ybtm    equal        ${ydimmin}-15

region  whole  block  -0.002  ${xdim} ${ybtm} ${ytop} 0.000000 ${zdim} units box
create_box  4  whole

region upper block INF INF -0.02 ${ydimmax} INF INF units box  
lattice   custom   ${alat} &
a1  %f    %f 	%f         &
a2  %f    %f    %f         &
a3  %f    %f    %f         &
%s
create_atoms 1 region upper

region lower block INF INF ${ydimmin} 0.02 INF INF units box
lattice   custom   ${alat} &
a1  %f    %f 	%f     &
a2  %f    %f    %f     &
a3  %f    %f    %f     &
%s
create_atoms 2 region lower

region uplayer  block INF INF ${yfreedimmax} ${ydimmax} INF INF units box
region lowlayer block INF INF ${ydimmin} ${yfreedimmin} INF INF units box

group upper region upper
group lower region lower

group uplayer  region uplayer
group lowlayer region lowlayer

group rigbody union uplayer lowlayer
group gbregion subtract all rigbody

set group lowlayer type 3
set group uplayer  type 4

pair_style   %s

# MEAM Kim
# pair_coeff  * *  ../lib_MgNdPb.meam  Mg  ../MgNd_para.meam  Mg  Mg  Mg  Mg

# EAM/alloy Poco
pair_coeff  * *  ../%s   Mg  Mg  Mg  Mg

neighbor 	     2.0      bin
displace_atoms upper    move   %f  0.0  %f  units box
delete_atoms   overlap  0.3  lower upper    # default

compute eng all pe/atom
compute eatoms gbregion reduce sum c_eng

fix f1 lowlayer setforce 0.0 0.0 0.0
fix f2 uplayer  aveforce 0.0 0.0 0.0

thermo  100
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
dump	 1 all custom 100 out/rel.chkpt.*  id type mass x y z c_eng

min_style cg
minimize 1e-15 1e-15 10000 10000

variable Ecoh equal   %f
variable Ntot equal   "count(gbregion)"
variable GBarea equal "lx*lz"
variable Ebulk equal  "v_Ecoh*count(gbregion)"
variable Etot equal   "c_eatoms"
variable Egb equal    "(c_eatoms-v_Ecoh*count(gbregion))/lx/lz*16.021"

print "Energy per atom =  ${Ecoh}"
print "Bulk energy ${Ebulk}"
print "Etot ${Etot}"
print "AreaXY = ${GBarea}"
print "N total atoms = ${Ntot}"
print "Egb = ${Egb}" #append ${f} screen yes
print "${Egb}" file lmp.dat screen yes
"""


inmdhcp = """
units metal
atom_style atomic
boundary  p  p  p 

variable  alat  equal    %f
variable  clat  equal    %f
variable  ysize equal    15 
variable  zsize equal    4 
variable  xdim  equal    %f
variable  zdim  equal    %f*${zsize}
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

variable ydimmin equal       -%d
variable ydimmax equal        %d
variable yfreedimmin equal   -%d
variable yfreedimmax equal    %d
variable extx    equal        0
variable ytop    equal        ${ydimmax}+30
variable ybtm    equal        ${ydimmin}-30
variable xtop    equal        ${xdim}+${extx}
variable xbtm    equal       -${extx}

region  whole  block  ${xbtm} ${xtop} ${ybtm} ${ytop} 0.000000 ${zdim} units box
create_box  4  whole

region upper block -0.01 ${xdim} -0.02 ${ydimmax} INF INF units box  
lattice   custom   ${alat} &
a1  %f    %f    %f         &
a2  %f    %f    %f         &
a3  %f    %f    %f         &
%s
create_atoms 1 region upper

region lower block -0.01 ${xdim} ${ydimmin} 0.02 INF INF units box
lattice   custom   ${alat} &
a1  %f    %f    %f     &
a2  %f    %f    %f     &
a3  %f    %f    %f     &
%s
create_atoms 2 region lower

region uplayer  block INF INF ${yfreedimmax} ${ydimmax} INF INF units box
region lowlayer block INF INF ${ydimmin} ${yfreedimmin} INF INF units box

group upper region upper
group lower region lower

group uplayer  region uplayer
group lowlayer region lowlayer

group rigbody union uplayer lowlayer
group gbregion subtract all rigbody

set group lowlayer type 3
set group uplayer  type 4

pair_style      %s
pair_coeff      * *  ../%s   Mg  Mg  Mg  Mg

neighbor        2.0  bin
displace_atoms  upper move  %f  0.0  %f  units box
delete_atoms    overlap  0.3 lower upper       # default

compute eng all pe/atom
compute eatoms  gbregion reduce sum c_eng

fix f1 lowlayer setforce 0.0 0.0 0.0
#fix f2 uplayer  aveforce 0.0 0.0 0.0

variable Ecoh equal %f
variable Ntot equal "count(gbregion)"
variable GBarea equal "lx*lz"
variable Ebulk equal "v_Ecoh*count(gbregion)"
variable Etot equal  "c_eatoms"
variable Egb equal   "(c_eatoms-v_Ecoh*count(gbregion))/lx/lz*16.021"

thermo    1000
thermo_style custom step temp etotal pxx pyy pzz v_Egb 
dump      1 all custom   5000  out/rel.chkpt.*  id type mass x y z c_eng 
variable  mtemp   equal  550 
fix       1       all    npt   temp   0.01   ${mtemp}   100  x  0  0  500  y  0  0  500 z 0 0 500  couple none  drag  1.0
restart   5000    restart/rst.*   
run       1000000

print "Energy per atom =  ${Ecoh}"
print "Bulk energy ${Ebulk}"
print "Etot ${Etot}"
print "AreaXY = ${GBarea}"
print "N total atoms = ${Ntot}"
print "Egb = ${Egb}" #append ${f} screen yes
print "${Egb}" file lmp.dat screen yes
"""
