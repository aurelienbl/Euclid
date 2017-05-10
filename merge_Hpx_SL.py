import os
import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt
import pymangle as pym
import healpy as hp


nside_out=4096
nside_tiny=4096*4

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
nside_L=1

# Je veux etre sur de prendre tous les tiles dans ce gros pixels, donc je cherche tous les fils a ndie=128
#m=nm.zeros(12*nside_L**2)
#m[L]=1
#hp.mollview(m, nest=True)




## la va etre le debut de la focntion qui prend en entree L. 
# verify that L is in dlarge

## something like if l is in dlarge then continu.


ng_gen=nm.log(nside_tiny/nside_out)/nm.log(2)

Npix_per_out=4**  ( nm.log(nside_tiny/nside_out)/nm.log(2))  





#nside_out=4
#nside_tiny=8

suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'

suff2='mask_V2'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'

suff3='mask_star_VIS_v2'
path_star='/data63/benoitl/Euclid/full_sky/'+suff3+'/'




coord='EQU'

map_out=nm.zeros(12*nside_out**2)+hp.UNSEEN
map_frac2=nm.zeros(12*nside_out**2)+hp.UNSEEN     # fraction of the subpixels that are seen at least 2 
map_frac3=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 3
map_frac4=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 4
map_tiny=nm.zeros(12*nside_tiny**2, dtype=nm.float128) +hp.UNSEEN


Nside_Super_large=1


for SL in range(12):
    print SL
    pixo=m_u.get_sons2([SL], 1,nside_out )
    pixt=m_u.get_sons2([SL], 1,nside_tiny)

    mo=hp.read_map(path_mask+'VIS_stars_Nobs_mean_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord, SL), nest=True)
    m2=hp.read_map(path_mask+'VIS_stars_Nobs_frac2_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord,  SL), nest=True)
    m3=hp.read_map(path_mask+'VIS_stars_Nobs_frac3_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord,  SL), nest=True)
    m4=hp.read_map(path_mask+'VIS_stars_Nobs_frac4_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord,  SL), nest=True)
    mt=hp.read_map(path_mask+'VIS_stars_Nobs_tiny_t.%d_%s_NESTED_SL%d.fits'%( nside_tiny,coord,  SL), nest=True)


    map_out[pixo]=mo[pixo]
    map_frac2[pixo]=m2[pixo]
    map_frac3[pixo]=m3[pixo]
    map_frac4[pixo]=m4[pixo]
    map_tiny[pixt]=mt[pixt]

    del(mo, m2,m3,m4,mt)
    


hp.write_map(path_mask+'VIS_stars_Nobs_mean_o.%d_t.%d_%s_NESTED.fits'%(nside_out, nside_tiny, coord),map_out, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac2_o.%d_t.%d_%s_NESTED.fits'%(nside_out, nside_tiny, coord),map_frac2, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac3_o.%d_t.%d_%s_NESTED.fits'%(nside_out, nside_tiny, coord),map_frac3, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac4_o.%d_t.%d_%s_NESTED.fits'%(nside_out, nside_tiny, coord),map_frac4, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_tiny_t.%d_%s_NESTED.fits'%( nside_tiny, coord),map_tiny, nest=True )




