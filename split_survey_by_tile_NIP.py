import os


import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt

import pymangle as pym


## Aims at generating a tiling for Euclid to create the Full Sky Mangle masks.  
## A setup is characterized by ra_min, ra_max, dec_min, dec_max, N_ra, N_dec, + a Mangle pixelisation level
## and a path.





suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'







suff2='mask_V2_NIP'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'
if not(os.path.isdir(path_mask)):
    os.mkdir(path_mask)




f=open(path_tiling+'params_tiling.json','r')
par=json.load( f)
f.close()


kmax=par["N_ra"]*par["N_dec"]

if not(os.path.isdir(path_mask+'mangle/')):
    os.mkdir(path_mask+'mangle/')




## path+'tiles/tiles_%d_%d.pol'%(kmax, i)))


VISw=pd.DataFrame.from_csv('/data63/benoitl/Euclid/Data_from_Joao/4400-nisp-footprint_WIDE_SURVEY.dat')

ra0=nm.array(VISw['RA0'])
ra1=nm.array(VISw['RA1'])
ra2=nm.array(VISw['RA2'])
ra3=nm.array(VISw['RA3'])

dec0=nm.array(VISw['DEC0'])
dec1=nm.array(VISw['DEC1'])
dec2=nm.array(VISw['DEC2'])
dec3=nm.array(VISw['DEC3'])



for i in range(1500, 1900):
    print  "doing tile ", i
    if not(os.path.isdir(path_mask+'mangle/tile_%d'%i)):
        os.mkdir(path_mask+'mangle/tile_%d'%i)


    ##load tile
    M=pym.Mangle(path_tiling+'tiles/tiles_%d_%d.pol'%(kmax, i))



    ## Get all the images that have at least one corner in that tile.
    ind=nm.where(  (M.contains(VISw['RA0'], VISw['DEC0'])==1) +(M.contains(VISw['RA1'], VISw['DEC1'])==1) + (M.contains(VISw['RA2'], VISw['DEC2'])==1)+ (M.contains(VISw['RA3'], VISw['DEC3'])==1))[0]
    


    N=ind.shape[0]
    if N!=0:
        f=open(path_mask+'mangle/tile_%d/subsample.edge'%i, 'w')
        
        for ii in range(N):
            s=mu.c2e2(ra0[ind[ii]], dec0[ind[ii]],  ra1[ind[ii]],dec1[ind[ii]], ra2[ind[ii]], dec2[ind[ii]], ra3 [ind[ii]], dec3[ind[ii]]  )
            print >>f, "edges %d"%ii
            print >>f, s[1]
        
        f.close()
 


        os.system('poly2poly -ie2 -k,1e-2  -m1e-8 %s  %s   '%(path_mask+'mangle/tile_%d/subsample.edge'%i, path_mask+'mangle/tile_%d/subsample.pol'%i)  )


        M=pym.Mangle( path_mask+'mangle/tile_%d/subsample.pol'%i)  
        #plt.figure()
        #plt.hist(M.areas, bins=100)



        os.system('pixelize -Ps0,%d -m1e-8  %s %s'%(par["pix"], path_mask+'mangle/tile_%d/subsample.pol'%i,  path_mask+'mangle/tile_%d/subsample_pixed.pol'%i ))
        os.system('snap  -m1e-8 -a0.000027 -b0.000027 -t0.000027  %s %s ' %( path_mask+'mangle/tile_%d/subsample_pixed.pol'%i, path_mask+'mangle/tile_%d/subsample_pixed_snapped.pol'%i ) )


        os.system('rasterize -T  -m1e-8 %s %s %s  '%( path_tiling+'tiles/tiles_%d_%d.pol'%(kmax, i),path_mask+'mangle/tile_%d/subsample_pixed_snapped.pol'%i ,  path_mask+'mangle/tile_%d/subsample_pixed_rasted.pol'%i      )    )
        
        
        os.system('balkanize -Ba  -m1e-8 %s %s ' %(  path_mask+'mangle/tile_%d/subsample_pixed_rasted.pol'%i,  path_mask+'mangle/tile_%d/subsample_pixed_balked.pol'%i ) )
        
        
