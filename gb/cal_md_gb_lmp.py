#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-11-01 01:54:26
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-06 15:34:53


class md_gb_lmp(object):

    def build_hcp_gb_lmp(self, ang=15, xdim=10, thk=(40, 50), disp=(0, 0, 0)):
        fid = open("in.gb", 'w')
        # pmat = self.hcp_til_mtx_z0001(ang)
        # nmat = self.hcp_til_mtx_z0001(-ang)
        # pmat = self.hcp_til_mtx_z1100(ang)
        # nmat = self.hcp_til_mtx_z1100(-ang)
        pmat = self.hcp_til_mtx_z1120(ang)
        nmat = self.hcp_til_mtx_z1120(-ang)
        fid.write(inhcp % (self.pot['ahcp'], self.pot['chcp'],
                           xdim,
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
                           self.pot['file'], disp[0], self.pot['ehcp']))

inhcp = """
units metal
atom_style atomic

variable  alat  equal    %f
variable  clat  equal    %f
variable  ysize equal    15 
variable  zsize equal    2
variable  xdim  equal    %f
variable  ydim  equal    sqrt(3)*${alat}*${ysize}
variable  zdim  equal    ${clat}*${zsize}
variable  sqa   equal    sqrt(3)
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

variable ydimmin equal    	 -%d
variable ydimmax equal   	    %d
variable yfreedimmin equal   -%d
variable yfreedimmax equal    %d

region  whole  block  -0.002  ${xdim} ${ydimmin} ${ydimmax} 0.000000 ${zdim} units box

create_box  4  whole

region upper block INF INF -0.02 ${ydimmax} INF INF units box  
lattice   custom   ${alat} &
a1  %f  	%f 		%f         &
a2  %f    %f    %f         &
a3  %f    %f    %f         &
%s
create_atoms 1 region upper

region lower block INF INF ${ydimmin} 0.02 INF INF units box
lattice   custom   ${alat} &
a1  %f  	%f 		%f     &
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

# ---------- Define Interatomic Potential ---------------------
pair_style   %s

# MEAM Kim
# pair_coeff  * *  ../lib_MgNdPb.meam  Mg  ../MgNd_para.meam  Mg  Mg  Mg  Mg

# EAM/alloy Poco
pair_coeff  * *  ../%s   Mg  Mg  Mg  Mg

neighbor 	 2.0  bin
neigh_modify delay 5 check yes

displace_atoms upper move  %f  0.0  0.0  units box

delete_atoms  overlap 0.3 lower upper    # default
#delete_atoms overlap 0.1 lower upper

compute eng all pe/atom
compute eatoms gbregion reduce sum c_eng

fix f1 lowlayer setforce 0.0 0.0 0.0
fix f2 uplayer  aveforce 0.0 0.0 0.0

thermo  100
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms
dump	 1 all custom 100 out/rel.chkpt.*  id type mass x y z c_eng
#restart 1 out/rel.restart

min_style cg
minimize 1e-12 1e-12 10000 10000

variable Ecoh equal %f
variable Ntot equal "count(gbregion)"
variable GBarea equal "lx*lz"
variable Ebulk equal "v_Ecoh*count(gbregion)"
variable Etot equal  "c_eatoms"
variable Egb equal   "(c_eatoms-v_Ecoh*count(gbregion))/lx/lz*16.021"

print "Energy per atom =  ${Ecoh}"
print "Bulk energy ${Ebulk}"
print "Etot ${Etot}"
print "AreaXY = ${GBarea}"
print "N total atoms = ${Ntot}"
print "Egb = ${Egb}" #append ${f} screen yes
print "${Egb}" file lmp.dat screen yes
"""
