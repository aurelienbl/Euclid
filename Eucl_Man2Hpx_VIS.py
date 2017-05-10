#! /usr/bin/env python



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


L=nm.int32(sys.argv[1])-1

print L


nside_out=4096
nside_tiny=8192*2

nside_ini=256

nside_L=1

coord='GAL'












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




FP=pym.Mangle(footprint_file)

ng_gen=nm.log(nside_tiny/nside_out)/nm.log(2)
Npix_per_out=4**  ( nm.log(nside_tiny/nside_out)/nm.log(2))  





map_out=nm.zeros(12*nside_out**2)+hp.UNSEEN
map_frac2=nm.zeros(12*nside_out**2)+hp.UNSEEN     # fraction of the subpixels that are seen at least 2 
map_frac3=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 3
map_frac4=nm.zeros(12*nside_out**2)+hp.UNSEEN    # fraction of the subpixels that are seen at least 4
map_tiny=nm.zeros(12*nside_tiny**2, dtype=nm.float128) +hp.UNSEEN



footprint=FP

map_out,map_frac2, map_frac3, map_frac4, map_tiny=em2h.perform_pixel_L(L, nside_L,nside_out, nside_tiny,nside_ini, coord, path_mask,path_star,footprint, map_out,map_frac2, map_frac3, map_frac4, map_tiny)

#map_out,map_frac2, map_frac3, map_frac4=em2h.perform_pixel_L(L, nside_L,nside_out, nside_tiny,nside_ini, coord, path_mask,path_star,footprint, map_out,map_frac2, map_frac3, map_frac4)

   ## map_out,map_frac2, map_frac3, map_frac4, map_tiny =perform_pixel_L(i, map_out,map_frac2, map_frac3, map_frac4, map_tiny) 
    
  

hp.write_map(path_mask+'VIS_stars_Nobs_mean_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord, L),map_out, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac2_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord, L),map_frac2, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac3_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord, L),map_frac3, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_frac4_o.%d_t.%d_%s_NESTED_SL%d.fits'%(nside_out, nside_tiny,coord, L),map_frac4, nest=True )
hp.write_map(path_mask+'VIS_stars_Nobs_tiny_t.%d_%s_NESTED_SL%d.fits'%( nside_tiny,coord, L),map_tiny, nest=True )




