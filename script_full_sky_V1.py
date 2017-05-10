
### OBSOLETE####



import numpy as nm
import matplotlib.pyplot as plt
import healpy as hp

import pandas as pd
import pymangle as pym
import os

from mangle import graphmask
import mangle_utils as mu

import os
from abl_mangle import mangle_plot as mp
from Man2Hpx import m_utils as m_u

path='/Users/benoitl/Documents/Post_doc/Euclid/Scripts/full_sky/'

########### first thing: split the sphere into rectangle !


dec_max= 89
dec_min=-89
ra_min=0
ra_max=360

N_ra=50
N_dec=50

k=0
for i in range(N_ra):
    for j in range(N_dec):
        k+=1
        
        ra0=ra_min+1.0* (ra_max-ra_min)/ N_ra *(i)
        ra1=ra_min+1.0* (ra_max-ra_min)/ N_ra *(i+1)

        
        dec0=dec_min+1.0* (dec_max-dec_min)/ N_dec *(j)
        dec1=dec_min+1.0* (dec_max-dec_min)/ N_dec *(j+1)

        
        f=open(path+'rects/try_%d.rect'%k, 'w')
        s=mu.c2r(ra0, ra1, dec0, dec1)
        print >>f, "rectangle %d"%k
        print >>f, s[1]
        f.close()
        os.system('poly2poly -ir1  -m1e-8 -q  %s %s'%(path+'rects/try_%d.rect'%k, path+'pols/try_%d.pol'%k))
        
kmax=k
## merge all the pols. 
list1=''
for i  in range(1, kmax+1):
    list1+=(path+'pols/try_%d.pol'%i)+' '

os.system('poly2poly -m1e-8  %s %s '%(list1, 'full_sphere_%d.pol'%kmax))
os.system('snap -m1e-8 -a0.0001 -b0.0001 -t0.0001  %s %s '%('full_sphere_%d.pol'%kmax, 'full_sphere_snapped_%d.pol'%kmax))
os.system('pixelize -m1e-8 -Ps0,8  %s %s '%('full_sphere_snapped_%d.pol'%kmax, 'full_sphere_snapped_pixelized_%d.pol'%kmax))

#### extracting the tiles
for i in range(1, 10):
    os.system('poly2poly -J%d,%d full_sphere_snapped_pixelized_2500.pol %s'%(i,i,path+'tiles/tiles_%d_%d.pol'%(kmax, i)))











if 1==2:

    f=open('try1.rect', 'w')
    s=mu.c2r(30, 40, 86, 87)
    i=5
    print >>f, "rectangle %d"%i
    print >>f, s[1]

    f.close()
    os.system('poly2poly -ir1  -m1e-8  try1.rect try1.pol2')


    f=open('try1.vert', 'w')
    s=mu.c2v(30, 86,40, 86, 40,87, 30,87)
    i=5
    print >>f, "vertices %d"%i
    print >>f, s[1]

    f.close()
    os.system('poly2poly -iv4  -m1e-8  try1.vert try1.vert.pol')





    VIS_w=pd.DataFrame.from_csv('4400-nisp-footprint_WIDE_SURVEY.dat')





    vw=VIS_w[ (VIS_w['RA0'] > 00) * (VIS_w['RA0'] < 60)    *( VIS_w['DEC0'] < 89 )   *(VIS_w['DEC0'] > -89)]  
    N=vw.shape[0]
    N


    f=open('subsample.edge', 'w')


    for i in range(N):
        #l=A[ind[0][i],:]
        s=c2e(vw['RA0'][vw.index[i]], vw['DEC0'][vw.index[i]] ,  vw['RA1'][vw.index[i]], vw['DEC1'][vw.index[i]]   ,  vw['RA2'][vw.index[i]], vw['DEC2'][vw.index[i]]   ,  vw['RA3'][vw.index[i]], vw['DEC3'][vw.index[i]]         )
        #s= c2e(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7])
        print >>f, "edges %d"%i
        print >>f, s[1]

    f.close()



    os.system('poly2poly -ie2 -k,1e-5  -m1e-8 subsample.edge  subsample.pol')


    M=pym.Mangle('subsample.pol')


    plt.hist(M.areas, bins=100)



    os.system('pixelize -Ps0,10 -m1e-8   subsample.pol subsample_pix.pol')
    os.system('snap  -m1e-8 -a0.000027 -b0.000027 -t0.000027  subsample_pix.pol subsample_snap.pol')
    os.system('balkanize -Ba  -m1e-8   subsample_snap.pol subsample_bal.pol')

    plt.figure()

    graphmask.plot_mangle_map('./subsample_bal.pol',linewidth=0, cmap=cm.jet, alpha=1, autoscale=True, bgcolor='grey')
    plt.savefig('./ex2_mangle.png')

    draw_colorbar


    m2048=mu.make_mask_EQU('./subsample_bal.pol', 2048)
    hp.gnomview(m2048,rot=[59.0, -31,0], nest=True, reso=.12, xsize=1200)
    plt.savefig('./ex1_hpx2048.png')
    m4096=mu.make_mask_EQU('./subsample_bal.pol', 4096)
    hp.gnomview(m4096,rot=[59.0, -31,0], nest=True, reso=.12, xsize=1200)
    plt.savefig('./ex1_hpx4096.png')




    M=pym.Mangle('subsample_bal.pol')
    plt.figure()
    plt.hist(M.weights, bins=20, weights=M.areas, cumulative=-1)
    plt.xlabel('#Observations', fontsize=18)
    plt.ylabel('Area (sq.deg)', fontsize=18)
    plt.savefig('./ex1_hist.png')
