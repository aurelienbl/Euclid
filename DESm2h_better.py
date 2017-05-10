#### This is DES stuff. Is just there so that I can take part of it and port them to Euclid.

import matplotlib.pyplot as plt
import healpy as hp
import numpy as nm
import pymangle as pym
import os
from Man2Hpx import m_utils as mu




def DESm2h(tolyfile,band, nside_out, nside_tiny, coord='GAL',weights=['maglims'], nside_large=8, nside_ini=128, project='OPS', tilefile='' ):
    """
    Routine designed to convert DES mangle mask to healpix masks without having the full combined masks
    The only combined mask that is required is the "tolygon polygons" file, that represent the footprint  of the mask at the tile level
    
    Inputs:
    
    tolyfile : toylogon file
    
    Parameters:
    band: band in g r i z Y
    nside_out: Nside of the final Healpix mask that are produced  (typically nside_out=2048 or 4096)
    nside_tiny: Nside of the finest resolution (typically nside_tiny=16384 or 32768)
    weights: list os strings that  indicate what maps are to be produced. Typically weights=['maglims', 'bitmask', 'time', 'weight']
    nside_large: very small nside (coarse resultion to minimize time. It is faster to run polyid on a big array of points than multiple polyid calls with smaller arrays ). Typically nside_large=8 or 16
    nside_ini : smallest nside needed tobe sure to probe all the tolygon file.  One tile is .5 sq deg. that requires probing the nside=128 pixels 
    coord: coordinante system. Should be only GAL or EQU
    project: Either ACT or OPS. used to find the mangle files on deslogin

    Outputs:

    N+1 healpix maps at resolution Nside_out, in Nested format and in the desired coordinates, where N is the number of weights considered. The first map (ind =0), is the detection fraction, i.e. the fraction of tiny pixel that are not 0 inside a final pixel.
    
   
    
    """
    if coord!='GAL' and coord!='EQU':
        print  "Please indidicate a correct coordinate system. Use either GAL or EQU"
        return -1

    nmaps=len(weights)
    maps=nm.zeros((12*nside_out**2, nmaps+1))+hp.UNSEEN


    ####################################################################################################
    ### 1. get at least one healpix pixel center within each tile of the mask. One tile is .5 sq deg. that requires probing the nside=128 pixels 
    ######################################################


    ## load ra, dec of healpix center for nside=nside_ini and galactic coordinates
    print 'Loading healpix center coordinates for nside= %d'%nside_ini
    
    if coord=='GAL':
        
        ra,dec= mu.get_hpx_coord_GAL_inradec(nside_ini)
    if coord=='EQU':
        ra,dec= mu.get_hpx_coord_radec(nside_ini)
    
    #Load tolygon files 
    print 'loading tolygons file'
    toly=pym.Mangle(tolyfile)
    print 'done'

    ##### get the weight of the nside=128 pixels against the tolygon mask. Tell whether pixels are in the mask
    w=toly.weight(ra,dec)


    ind=nm.where(w !=0)[0]


    ### ind contains the pixels that are in the mask. to make sure that we are indeed having all the pixels that have even the tiniest overlap with the mask, we take twice the neighbouring pixels

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
    ### 2. consider even larger pixels to minimize io time. it is faster to run polyid on a big array of points than multiple polyid calls with smaller arrays
    ##########################################################################################
    dlarge=nm.int32(d/4**(nm.log(nside_ini/nside_large)/nm.log(2)))
    dlarge=nm.unique(dlarge)

    print "There are %d large nside=%d pixels in this mask'"%(dlarge.shape[0], nside_large)
    # dlarge now contains the ids of the nside= nside_large pixels

    ############################################################################################
    ### 3. loading the table with the mangle_run, coadd_run ad coaddtile_id
    ############################################################################################

    #### Need to add something to get the tiles.dat files


#    ff=open('/home/benoitl/dev2/scripts/mangle_Y1T18/Y1A1_S82_tiles.dat')
    ff=open('/home/benoitl/y1a1_spt_tiles.dat')
  ##  ff=open(tilefile)
    lines=ff.readlines()
    N=len(lines)-1
    dum1=[]
    for i in range(1,N+1):
        dum=lines[i].strip().split(',')
        dum1.append(dum)
    tile_tab=nm.array(dum1)
    ff.close()

    #### tile_tab contains all the info af the tiles in the mask: coaddtile_id, tilename, mangle_run, coadd_run. Those should be needed to retreive the individual tiles masks



    npix_per_fat=nm.int(4**(nm.log(nside_tiny/nside_out)/nm.log(2)))  ## number of tiny pixels per pixels of the final healpix map.

    kk=0
    for lpixid in dlarge:    ##### loop on the large pixels
        kk+=1
        print '%d / %d big pixels'%(kk, dlarge.shape[0])

        ## get all the pixels of the final map for the large pixel lpixid 
        fat=mu.get_sons2(nm.array([lpixid]), nside_large,nside_out)

        ##these are the ipix of all the tiny hpx pixels within the fat 
        sons=mu.get_sons2(fat, nside_out, nside_tiny)
        jmaps=nm.zeros((len(sons), nmaps+1))+hp.UNSEEN 


        ####################################################################
        ### Need to transform back the coord of the sons into EQU to be able to match them agains the mangle mask

        Npixel=sons.shape[0]   ### Npixel is the number of tiny pixels with a large one.

        R=hp.Rotator(coord='gc')

        theta,phi=hp.pix2ang(nside_tiny, sons, nest=True)
        theta=nm.pi/2 -theta
        if coord=='GAL':
            ra, dec= R(phi*180/nm.pi, theta*180/nm.pi, lonlat=True)
            ra=nm.mod(ra, 360)
        if coord=='EQU':
            ra=phi*180/nm.pi
            dec=theta*180/nm.pi

            ra=nm.mod(ra, 360)
            ra,dec=mu.perturb_radec(nside_tiny, ra,dec)



        #### Need to split the sons pixels according to which tile thy fall in.
        p=toly.polyid(ra, dec)
        a,b,c=nm.unique(p, return_index=True, return_inverse=True)

        print 'There are %d tiles in this large pixel of the sky in which we will look one by one'%(a.shape[0])
        
        #####For each tile within this large pixel lpixid interoggate the mangle masks at the position of the tiny pixels.

        for i in  range(a.shape[0]):
            t1=a[i]
            if t1!=-1: 
                ind_1= nm.where(tile_tab[:,0]=='%s'%t1)[0]  ## get the line in tile_tab that correspond to that tile
            
                mangle_run=tile_tab[ind_1,2][0]
                coadd_run=tile_tab[ind_1,3][0]
                del(ind_1)
                # t1 is the tileid of the firsttile in which there are some pixles.
                #ind1=where(c==i) are the indices of those pixels
                indd=nm.where(c==i) 
                verbose=False

                # loading the mangle mask for that first tile. The path of the files with have to be modified when ported to the deslogin
                if verbose:
                    print 'loading file /archive_data/Archive/%s/mangle/%s/mangle/%s_holymolys_weight_%s.pol'%(project, mangle_run, coadd_run,band)

                if os.path.isfile('/archive_data/Archive/%s/mangle/%s/mangle/%s_holymolys_weight_%s.pol'%(project,mangle_run, coadd_run,band))==False:
                    print '/archive_data/Archive/%s/mangle/%s/mangle/%s_holymolys_weight_%s.pol  not found. Skipping'%(project,mangle_run, coadd_run,band)
                else:

                    nmg1=pym.Mangle('/archive_data/Archive/%s/mangle/%s/mangle/%s_holymolys_weight_%s.pol'%(project,mangle_run, coadd_run,band))

                    ##reweight these polygons with 1 2 3 4  Npolys, so that we can select the correct values in the weight files.
                    nm.savetxt('range_npolys%d_%s'%(nside_tiny,band), nm.arange(1, 1+nmg1.get_npoly()))
                    nmg1.read_weights('range_npolys%d_%s'%(nside_tiny, band))
                    os.system('rm range_npolys%d_%s'%(nside_tiny, band))
                    ids=nm.int32(nmg1.weight(ra[indd], dec[indd]))
                    ### ids -1 would be the index of the polygon in the weight files. There need to be a 1 shift as weight=0 means that the ra,dec is not in the mask.

                    we=nm.loadtxt('/archive_data/Archive/%s/mangle/%s/mangle/%s_molys_%s.bitmask'%(project,mangle_run, coadd_run,band))



                    we[we>0]=-2
                    we[we==0]=1
                    we[we==-2]=0 ### These thres lines is to set to 1 the good regions and 0 the regions with bitmask >0


                    w1=ra[indd].copy()+hp.UNSEEN

                    w1=we[ids-1]
                    w1[ids==0]=hp.UNSEEN

                


    ##                w1= nmg1.weight(ra[indd], dec[indd])
                    jmaps[indd,0]=nm.float64(w1)

                    #                jweight[indd]=nm.float64(w1)
                    del(w1,we)
                    it=0
                    for  wei in weights:
                        it+=1

                        weiarr=nm.loadtxt('/archive_data/Archive/%s/mangle/%s/mangle/%s_molys_%s.%s'%(project, mangle_run, coadd_run,band, wei))

                        w3=weiarr[ids-1]
                        w3[ids==0]=hp.UNSEEN

                        jmaps[indd, it]=nm.float64(w3)

                        del(w3)
                        ###fin boucle sur tile 

                        ### load maglim
    #                    nmg1.read_weights('/archive_data/Archive/%s/mangle/%s/mangle/%s_molys_%s.%s'%(project,mangle_run, coadd_run, band,wei))
     #                   w1= nmg1.weight(ra[indd], dec[indd])
      #                  jmaps[indd, it]=nm.float64(w1)
                        #                    jmaglim[indd]=nm.float64(w1)
       #                 del(w1)
                        ###fin boucle sur tile 
                    del(ids)
        jjmaps=nm.reshape(jmaps, (Npixel/npix_per_fat, npix_per_fat,nmaps+1))
    
        mm=nm.mean(jjmaps, axis=1, dtype=nm.float64)
        

        mm2=mm.copy()*0.0+ hp.UNSEEN
        mfrac=mm2[:,0].copy()*0.0
        
        for j in range(1,nmaps+1):
            
            for i in range(Npixel/npix_per_fat):    ## tjis the number of fat pixels within  a large one
                ### a where to get the pixels not unseen
                idd=nm.where(jjmaps[i, :, j]!=hp.UNSEEN)[0]
                if len(idd)>0:
                    mm2[i, j]=nm.mean(jjmaps[i, idd, j])
                 

        for i in range(Npixel/npix_per_fat):
            if len(nm.where(jjmaps[i,:,0]==hp.UNSEEN)[0])==npix_per_fat:
                mfrac[i]=hp.UNSEEN
            else:
                mfrac[i]=nm.sum(jjmaps[i,:,0]>0)/(1.0*npix_per_fat)
       # mfrac2=nm.mean(jjmaps[:,:,0][jjmaps[:,:,0]>0], axis=1)


#        mfrac=nm.sum(jjmaps[:,:,0]>0, axis=1)/(1.0*npix_per_fat)
        maps[fat,1:]=mm2[:,1:]
        maps[fat,0]=mfrac
        del( sons, fat, jjmaps,jmaps,mm, mfrac,)
        del(p, a,b,c)
    return maps











#root_dir= '/Users/benoitl/Documents/Post_doc/DES/mangle_healpix/'

#man_dir='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/coadd/20130320000001_DES0219-0541//mangle/'
#man_file='20130320000001_DES0219-0541_holymolys_weight_r.pol'

## directory of the combined mangle directory where the tolys.pol file is. 
#man_dir='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82/'
#man_file='y1p1_s82_tolys.pol'


##tolyfile='/Users/benoitl/Documents/Post_doc/DES/DESdata/OPS/Final_masks/mangle/y1p1_s82/y1p1_s82_tolys.pol'

#tolyfile=man_dir+man_file

#nside_large=8 ## Healpix map is done in each of these big pixels
#nside_ini=128 ## Initial pixels udes to find where the mangle mask is, should be failty big, to have a dense ampling of the tolygon polygon file
#nside_map_final=1024  ## Final resolution of the output healpix maps
#nside_tiny=4096   ### finer scale on which the averaging is done



#weights=['maglims', 'bitmask', 'time', 'weight']


#maps_G=DESm2h_gal(tolyfile, nside_map_final, nside_tiny, weights=weights, nside_large=nside_large, nside_ini=128 )
#maps_C=m2h_equ(tolyfile, nside_map_final, nside_tiny, weights=weights, nside_large=16, nside_ini=128 )









