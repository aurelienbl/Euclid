
import os
import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt
import pymangle as pym
import healpy as hp
import Eucl_Man2Hpx as em2h 
import sys
import pyfits as pyf
import time



coord='GAL'









### mask stuff


suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'

suff2='mask_V2_NIP'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'

suff3='mask_star_J_v1'
path_star='/data63/benoitl/Euclid/full_sky/'+suff3+'/'



f=open(path_tiling+'params_tiling.json','r')
par=json.load( f)
f.close()

footprint_file=path_mask+'footprint.pol'
kmax=par["N_ra"]*par["N_dec"]




## catalog stuff

catalogf='/data63/benoitl/Euclid/MICE2/82.fits'

t0=time.time()


h=pyf.open(catalogf)
ram=h[1].data.field('ra_gal_mag')
decm=h[1].data.field('dec_gal_mag')
#dec=90.-dec
ram=ram+180


i_good=nm.where(nm.isfinite(ram))
i_bad=nm.where(nm.isfinite(ram)==False)


ram[i_bad]=0.0
decm[i_bad]=0.0

print "time to load the MICE2 catalog,and add 180 to RA: %f sec "%(time.time()-t0)

#ram=ram[i_good]
#decm=decm[i_good]


## changing coords as the mangle masks are in ECL.
t0=time.time()
R=hp.Rotator(coord='CE')
ram1, decm1= R(ram, decm, lonlat=True)
ram1=nm.mod(ram1, 360)

print "time to Rotate the MICE2 catalog from C to E: %f sec "%(time.time()-t0)



if 1==2:
    plt.figure(1)
    plt.hist2d(ram[::100], decm[::100], bins=50)



    ra1, dec1= R(ram[::100], decm[::100], lonlat=True)
    ra1=nm.mod(ra1, 360)
    
    
    plt.figure(2)
    plt.hist2d(ra1, dec1, bins=50)
    

pathm='/data63/benoitl/Euclid/MICE2/'

FP=pym.Mangle(footprint_file)


if 1==2:
    
    t0=time.time()
    p=FP.polyid(ram1[0:10000000], decm1[0:10000000])
    print ram1[0:10000000].shape
    print time.time()-t0
    


    t0=time.time()
    p=FP.polyid(ram1[0:100000000], decm1[0:100000000])
    print ram1[0:100000000].shape
    print time.time()-t0



    t0=time.time()
    p=FP.polyid(ram1[0:200000000], decm1[0:200000000])
    print ram1[0:200000000].shape
    print time.time()-t0
    



    t0=time.time()
    p=FP.polyid(ram1[0:300000000], decm1[0:300000000])
    print ram1[0:300000000].shape
    print time.time()-t0
    


t0=time.time()
p=FP.polyid(ram1, decm1)
print ram1.shape
print "time to run the catalog through the footprint: %f sec "%(time.time()-t0)


    
N_obs=nm.zeros(len(ram1), dtype=nm.float128)-1
star_vis=nm.zeros(len(ram1), dtype=nm.float128)


t0=time.time()

print 't3'
a,b,c=nm.unique(p, return_index=True, return_inverse=True)
print 't4'
print "time find the unique tiles: %f sec "%(time.time()-t0)


    ## Number of tiles in this large pixel
Ntile=a.shape[0]


t0=time.time()

print 'There are %d tiles in this large pixel of the sky in which we will look one by one'%(Ntile)

for i in range(a.shape[0]):
    id_tile=a[i]
    if id_tile!=-1:
        if nm.mod(i, 50)==0:
            print 'doing tile # %d'%i
        indd=nm.where(c==i)
        mangle_mask_file=path_mask+'mangle/tile_%d/subsample_pixed_balked.pol'%id_tile
        Mg=pym.Mangle(mangle_mask_file)
        DD=Mg.polyid_and_weight(ram1[indd],decm1[indd])
        outside=nm.where(DD[0]==-1)[0]

        N_obs[indd]=DD[1]

        star_mask_file=path_star+'mangle/tile_%d/subsamplecirc_pixed_rasted.pol'%(id_tile)
        if  os.path.isfile(star_mask_file):
        #    print 'there are stars here'
            Ms=pym.Mangle(star_mask_file)
            EE=Ms.polyid_and_weight(ram1[indd],decm1[indd])

            star_vis[indd]=EE[1]

print "time to run the catalog through the  whole mask: %f sec "%(time.time()-t0)


N_obs[i_bad]=-2
star_vis[i_bad]=-2

t0=time.time()
em2h.write_indices_and_values(h[1].data.field('unique_gal_id'),N_obs ,'/data63/benoitl/Euclid/MICE2/Nobs_82_J.fits')
em2h.write_indices_and_values(h[1].data.field('unique_gal_id'),star_vis ,'/data63/benoitl/Euclid/MICE2/star_flag_82_J.fits')
print "time to save the results : %f sec "%(time.time()-t0)
