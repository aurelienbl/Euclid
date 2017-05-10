import os
import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt
import pymangle as pym
import healpy as hp


nside_out=1024
nside_tiny=1024*2

suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'

suff2='mask_V2'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'

suff3='mask_star_VIS_v2'
path_star='/data63/benoitl/Euclid/full_sky/'+suff3+'/'


f=open(path_tiling+'params_tiling.json','r')
par=json.load( f)
f.close()

footprint_file=path_mask+'footprint.pol'
kmax=par["N_ra"]*par["N_dec"]



# Je veux faire un gros pixel L a nside =8, 16, 32 ...
nside_L=8

# Je veux etre sur de prendre tous les tiles dans ce gros pixels, donc je cherche tous les fils a ndie=128
#m=nm.zeros(12*nside_L**2)
#m[L]=1
#hp.mollview(m, nest=True)



#A=m_u.get_sons2([L], 8, 128)
#r,d=hp.pix2ang(128, A, nest=True)
nside_ini=256

print 'Loading healpix center coordinates for nside= %d'%nside_ini
coord='EQU'

if coord=='GAL':
    ra,dec= m_u.get_hpx_coord_GAL_inradec(nside_ini)

if coord=='EQU':
    ra,dec= m_u.get_hpx_coord_radec(nside_ini)





    
#A me donne la liste des pixels a nside =128 dan le gros pixels L.

## checke que au moins 1 de ces pixels est dans le footprint.

FP=pym.Mangle(footprint_file)

B=FP.weight(ra,dec)
    #w=toly.weight(ra,dec)


ind=nm.where(B !=0)[0]


    ### ind contains the pixels that are in the mask. to make sure that we are indeed having all the pixels that have even the tiniest ove\
#rlap with the mask, we take twice the neighbouring pixels                                                                                  

a=hp.get_all_neighbours(nside_ini,ind, nest=True)
a=a[a>0]
b=nm.concatenate([a, ind])
b=nm.unique(b)
c=hp.get_all_neighbours(nside_ini,b, nest=True)
c=c[c>0]
c=nm.concatenate([c,b])
d=nm.unique(c)

    ### d now contains all the mid-scale pixels that we are going to consider.                                                             



    ###########################################################################################                                            
    ### 2. consider even larger pixels to minimize io time. it is faster to run polyid on a big array of points than multiple polyid calls\
# with smaller arrays                                                                                                                       
    ##########################################################################################                                             
dlarge=nm.int32(d/4**(nm.log(nside_ini/nside_L)/nm.log(2)))
dlarge=nm.unique(dlarge)




## la va etre le debut de la focntion qui prend en entree L. 
# verify that L is in dlarge

## something like if l is in dlarge then continu.


ng_gen=nm.log(nside_tiny/nside_out)/nm.log(2)

Npix_per_out=4**  ( nm.log(nside_tiny/nside_out)/nm.log(2))  




def perform_pixel_L(L, map_out,map_frac2, map_frac3, map_frac4, map_tiny):

    print  'Doing Large pixel %d at nside %s'%(L, nside_L)


    sons1=m_u.get_sons2(nm.array([L]), nside_L,nside_out)  ## these are the pixels in the output map (nside=2048, 4096) that are in the large pixel L
    sons2=m_u.get_sons2(sons1, nside_out, nside_tiny)  ## these are the smallest pixel that will subsample each output pixel

    sons2_values=nm.zeros((len(sons2),1))+hp.UNSEEN     ## this will be the values of the mangle mask at the centers of the tiny pixels

    Npixel=sons2.shape[0]   ### Npixel is the number of tiny pixels with a large one.

    R=hp.Rotator(coord='gc')

    theta,phi=hp.pix2ang(nside_tiny, sons2, nest=True)
    theta=nm.pi/2 -theta

    if coord=='GAL':
        ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
        ra=nm.mod(ra, 360)

    if coord=='EQU':
        ra=phi*180/nm.pi
        dec=theta*180/nm.pi
        ra=nm.mod(ra, 360)
        ra,dec=m_u.perturb_radec(nside_tiny, ra,dec)



    #### Need to split the sons2 pixels according to which tile thy fall in.

    p=FP.polyid(ra, dec)
    a,b,c=nm.unique(p, return_index=True, return_inverse=True)


    ## Number of tiles in this large pixel
    Ntile=a.shape[0]

    print 'There are %d tiles in this large pixel of the sky in which we will look one by one'%(Ntile)
    i=0
    for i in  range(a.shape[0]):
        id_tile=a[i]
        if id_tile!=-1:
       ## id of the tile under consideration
            print 'doing tile # %d'%i
            indd=nm.where(c==i)      ### That gives the tiny pixels that are in tile  i 
    ##if 1!=-1:    Ok this must be the case sometime.Let's try wihtou for now.
    ### Load the mangle mask in this tile
            mangle_mask_file=path_mask+'mangle/tile_%d/subsample_pixed_balked.pol'%id_tile
            Mg=pym.Mangle(mangle_mask_file)
            DD=Mg.polyid_and_weight(ra[indd],dec[indd])
            outside=nm.where(DD[0]==-1)[0]
            DD[1][outside]=hp.UNSEEN
                # check that there are some stars in this tile

            star_mask_file=path_star+'mangle/tile_%d/subsamplecirc_pixed_rasted.pol'%(id_tile)
            if  os.path.isfile(star_mask_file):
                print 'there are stars here'
                Ms=pym.Mangle(star_mask_file)
                EE=Ms.polyid_and_weight(ra[indd],dec[indd])
                
                value=DD[1]*(1- EE[1])
                value[outside]=hp.UNSEEN
                map_tiny[sons2[indd]]=value
                sons2_values[indd[0],0]=value

                #map_tiny[sons2[indd]]=DD[1]*(1- EE[1])
                #sons2_values[indd[0],0]=DD[1]*(1- EE[1])
            else:
                map_tiny[sons2[indd]]=DD[1]
                sons2_values[indd[0],0]=DD[1]



    for j in range(len(sons1)):
        ## subpixels
        subpix=sons2_values[j*4**ng_gen:(j+1)*4**ng_gen]
        nb_seen=len(nm.where(subpix!=hp.UNSEEN)[0])
        if nb_seen!= 0:  #Npix_per_out
            map_out[sons1[j]]=nm.mean(subpix[nm.where(subpix!=hp.UNSEEN)[0]])
            map_frac2[sons1[j]]=len (nm.where(subpix >=2  )[0])/Npix_per_out
            map_frac3[sons1[j]]=len (nm.where(subpix >=3  )[0])/Npix_per_out
            map_frac4[sons1[j]]=len (nm.where(subpix >=4  )[0])/Npix_per_out



    

    return map_out,map_frac2, map_frac3, map_frac4, map_tiny




map_out=nm.zeros(12*nside_out**2)+hp.UNSEEN
map_frac2=nm.zeros(12*nside_out**2)+hp.UNSEEN     # fraction of the subpixels that are seen at least 2 
map_frac3=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 3
map_frac4=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 4
map_tiny=nm.zeros(12*nside_tiny**2, dtype=nm.float128) +hp.UNSEEN






for i in range(0,12* nside_L**2):
    print  'i=%d/%d'%(i, 12*nside_L**2)
    map_out,map_frac2, map_frac3, map_frac4, map_tiny =perform_pixel_L(i, map_out,map_frac2, map_frac3, map_frac4, map_tiny) 
    
    if nm.mod(i,16)==0:
        hp.write_map(path_mask+'VIS_stars_Nobs_mean_o.%d_t.%d_ECL_NESTED_%d.fits'%(nside_out, nside_tiny, i),map_out, nest=True )
        hp.write_map(path_mask+'VIS_stars_Nobs_frac2_o.%d_t.%d_ECL_NESTED_%d.fits'%(nside_out, nside_tiny, i),map_frac2, nest=True )
        hp.write_map(path_mask+'VIS_stars_Nobs_frac3_o.%d_t.%d_ECL_NESTED_%d.fits'%(nside_out, nside_tiny, i),map_frac3, nest=True )
        hp.write_map(path_mask+'VIS_stars_Nobs_frac4_o.%d_t.%d_ECL_NESTED_%d.fits'%(nside_out, nside_tiny, i),map_frac4, nest=True )
        hp.write_map(path_mask+'VIS_stars_Nobs_tiny_t.%d_ECL_NESTED_%d.fits'%( nside_tiny, i),map_tiny, nest=True )

        

hp.write_map(path_mask+'VIS_stars_Nobs_mean_o.%d_t.%d_ECL_NESTED.fits'%(nside_out, nside_tiny),map_out, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac2_o.%d_t.%d_ECL_NESTED.fits'%(nside_out, nside_tiny),map_frac2, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac3_o.%d_t.%d_ECL_NESTED.fits'%(nside_out, nside_tiny),map_frac3, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac4_o.%d_t.%d_ECL_NESTED.fits'%(nside_out, nside_tiny),map_frac4, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_tiny_t.%d_ECL_NESTED.fits'%( nside_tiny),map_tiny, nest=True )




