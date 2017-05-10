import os
import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt
import pymangle as pym
import healpy as hp
import pyfits as pyf

def write_indices_and_values(ind, val, filename):
   
    
    table=[pyf.Column(name='INDICES',format='1K',array=ind),pyf.Column(name='VALUES',format='1D',array=val) ]
    #print table
    tbhdu=pyf.new_table(table)
    tbhdu.writeto(filename, clobber=True)

def read_indices_and_values(filename):
   

    h=pyf.open(filename)
    n_val=len(h[1].data.field(0))
    #indval=nm.zeros((2, n_val))
    #indval[0, :]=h[1].data.field(0)
    #indval[1, :]=h[1].data.field(1)
   
    return h[1].data.field(0), h[1].data.field(1)



def perform_pixel_L_V2(L, nside_L,nside_out, nside_tiny, nside_ini, coord, path_mask,path_star,footprint ):# map_out,map_frac2, map_frac3, map_frac4, map_tiny=''):
    ''' no need to provides initial maps. will return fits file with the list of pixid in nested and asociated values. Everything shoudl get MUCH faster!

    '''


    ng_gen=nm.log(nside_tiny/nside_out)/nm.log(2)
    #  the primary input is in ECL coordinate
    Npix_per_out=4**  ( nm.log(nside_tiny/nside_out)/nm.log(2))  
    print  'Doing Large pixel %d at nside %s'%(L, nside_L)
    sons1=m_u.get_sons2(nm.array([L]), nside_L,nside_out)  ## these are the pixels in the output map (nside=2048, 4096) that are in the large pixel L
    sons2=m_u.get_sons2(sons1, nside_out, nside_tiny)  ## these are the smallest pixel that will subsample each output pixel

    sons2_values=nm.zeros((len(sons2),1))+hp.UNSEEN     ## this will be the values of the mangle mask at the centers of the tiny pixels


    Npixel2=sons2.shape[0]   ### Npixel is the number of tiny pixels with a large one.
    Npixel1=sons1.shape[0] 

    map_out=nm.zeros((len(sons1),1))+hp.UNSEEN  
    map_frac2=nm.zeros((len(sons1),1))+hp.UNSEEN  
    map_frac3=nm.zeros((len(sons1),1))+hp.UNSEEN 
    map_frac4=nm.zeros((len(sons1),1))+hp.UNSEEN  
    #map_tiny
    
    print 't1'


    theta,phi=hp.pix2ang(nside_tiny, sons2, nest=True)
    theta=nm.pi/2 -theta

    if coord=='GAL':
        R=hp.Rotator(coord='ge')
        ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
        ra=nm.mod(ra, 360)

    if coord=='EQU':
        R=hp.Rotator(coord='ce')
        ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
        ra=nm.mod(ra, 360)


    if coord=='ECL':
        ra=phi*180/nm.pi
        dec=theta*180/nm.pi
        ra=nm.mod(ra, 360)
        ra,dec=m_u.perturb_radec(nside_tiny, ra,dec)

    print 't2'
    #### Need to split the sons2 pixels according to which tile thy fall in.

    p=footprint.polyid(ra, dec)
    print 't3'
    a,b,c=nm.unique(p, return_index=True, return_inverse=True)
    print 't4'
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
                #if map_tiny!='':
                 #   map_tiny[sons2[indd]]=value
                sons2_values[indd[0],0]=value

                #map_tiny[sons2[indd]]=DD[1]*(1- EE[1])
                #sons2_values[indd[0],0]=DD[1]*(1- EE[1])
            else:
                #if map_tiny!='':
                 #   map_tiny[sons2[indd]]=DD[1]
                sons2_values[indd[0],0]=DD[1]



    for j in range(len(sons1)):
        ## subpixels
        subpix=sons2_values[j*4**ng_gen:(j+1)*4**ng_gen]
        nb_seen=len(nm.where(subpix!=hp.UNSEEN)[0])
        if nb_seen!= 0:  #Npix_per_out
            map_out[j]=nm.mean(subpix[nm.where(subpix!=hp.UNSEEN)[0]])
            map_frac2[j]=len (nm.where(subpix >=2  )[0])/Npix_per_out
            map_frac3[j]=len (nm.where(subpix >=3  )[0])/Npix_per_out
            map_frac4[j]=len (nm.where(subpix >=4  )[0])/Npix_per_out


   # if map_tiny!='':           
    return sons1, sons2, map_out,map_frac2, map_frac3, map_frac4, sons2_values
    #else:
     #   return map_out,map_frac2, map_frac3, map_frac4




def perform_pixel_L(L, nside_L,nside_out, nside_tiny, nside_ini, coord, path_mask,path_star,footprint, map_out,map_frac2, map_frac3, map_frac4, map_tiny=''):

    ng_gen=nm.log(nside_tiny/nside_out)/nm.log(2)
    #  the primary input is in ECL coordinate
    Npix_per_out=4**  ( nm.log(nside_tiny/nside_out)/nm.log(2))  
    print  'Doing Large pixel %d at nside %s'%(L, nside_L)
    sons1=m_u.get_sons2(nm.array([L]), nside_L,nside_out)  ## these are the pixels in the output map (nside=2048, 4096) that are in the large pixel L
    sons2=m_u.get_sons2(sons1, nside_out, nside_tiny)  ## these are the smallest pixel that will subsample each output pixel

    sons2_values=nm.zeros((len(sons2),1))+hp.UNSEEN     ## this will be the values of the mangle mask at the centers of the tiny pixels
    Npixel=sons2.shape[0]   ### Npixel is the number of tiny pixels with a large one.


    theta,phi=hp.pix2ang(nside_tiny, sons2, nest=True)
    theta=nm.pi/2 -theta

    if coord=='GAL':
        R=hp.Rotator(coord='ge')
        ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
        ra=nm.mod(ra, 360)

    if coord=='EQU':
        R=hp.Rotator(coord='ce')
        ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
        ra=nm.mod(ra, 360)


    if coord=='ECL':
        ra=phi*180/nm.pi
        dec=theta*180/nm.pi
        ra=nm.mod(ra, 360)
        ra,dec=m_u.perturb_radec(nside_tiny, ra,dec)

    #### Need to split the sons2 pixels according to which tile thy fall in.

    p=footprint.polyid(ra, dec)
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
                if map_tiny!='':
                    map_tiny[sons2[indd]]=value
                sons2_values[indd[0],0]=value

                #map_tiny[sons2[indd]]=DD[1]*(1- EE[1])
                #sons2_values[indd[0],0]=DD[1]*(1- EE[1])
            else:
                if map_tiny!='':
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


    if map_tiny!='':           
        return map_out,map_frac2, map_frac3, map_frac4, map_tiny
    else:
        return map_out,map_frac2, map_frac3, map_frac4

