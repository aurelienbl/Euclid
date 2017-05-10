import pyfits as pyf
import os


import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt

import pymangle as pym
import healpy as hp


stars=pyf.open('/data63/benoitl/Euclid/2MASS_i_head.fits')[1].data



suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'



suff2='mask_star_J_v1'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'
if not(os.path.isdir(path_mask)):
    os.mkdir(path_mask)

f=open(path_tiling+'params_tiling.json','r')
par=json.load( f)
f.close()


kmax=par["N_ra"]*par["N_dec"]

if not(os.path.isdir(path_mask+'mangle/')):
    os.mkdir(path_mask+'mangle/')



##### Rotate those coordonates into Ecliptic
 

R=hp.Rotator(coord='GE', deg=True)
inds=nm.where((stars['i_AB']>4  ) * (stars['i_AB']<16  )   )[0]
Ns=inds.shape[0]; Ns
stars2=stars[inds]
rad=mu.radius_J(stars2['i_AB'])
rd0=R(stars2['GRA'], stars2['GDEC'], lonlat=True)
rd0[0]=nm.mod(rd0[0], 360)


 ## lets work with tile i
for i in range(2300, kmax+1):
    print  "doing tile ", i
    if not(os.path.isdir(path_mask+'mangle/tile_%d'%i)):
        os.mkdir(path_mask+'mangle/tile_%d'%i)


    ## get all the star center that are into tile i
    M=pym.Mangle(path_tiling+'tiles/tiles_%d_%d.pol'%(kmax, i))
    ind=nm.where(  M.contains(rd0[0], rd0[1])==1)[0]

    N=ind.shape[0]

    if N!=0:
        f=open(path_mask+'mangle/tile_%d/subsample.circ'%i, 'w')

        for ii in range(N):
            print >>f, rd0[0][ind[ii]],rd0[1][ind[ii]] , rad[ind[ii]]*0.1/60/60.


        f.close()
        os.system('poly2poly  -m1e-8 -ic   %s %s '%(path_mask+'mangle/tile_%d/subsample.circ'%i, path_mask+'mangle/tile_%d/subsamplecirc.pol'%i)  )
        os.system('pixelize -Ps0,%d -m1e-8  %s %s'%(par["pix"], path_mask+'mangle/tile_%d/subsamplecirc.pol'%i,  path_mask+'mangle/tile_%d/subsamplecirc_pixed.pol'%i ))
        os.system('snap  -m1e-8 -a0.000027 -b0.000027 -t0.000027  %s %s ' %( path_mask+'mangle/tile_%d/subsamplecirc_pixed.pol'%i, path_mask+'mangle/tile_%d/subsamplecirc_pixed_snapped.pol'%i ) )
        os.system('rasterize -T  -m1e-8 %s %s %s  '%( path_tiling+'tiles/tiles_%d_%d.pol'%(kmax, i),path_mask+'mangle/tile_%d/subsamplecirc_pixed_snapped.pol'%i ,  path_mask+'mangle/tile_%d/subsamplecirc_pixed_rasted.pol'%i      )    )



