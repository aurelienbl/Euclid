import os


import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u
import pandas as pd
import numpy as nm
import matplotlib.pyplot as plt
import pymangle as pym


## One a mangle mask is generated, this script aggregates the tiles which are not empty.




suff1='tiling_V2_pix10'
path_tiling='/data63/benoitl/Euclid/full_sky/'+suff1+'/'

suff2='mask_V2_NIP'
path_mask='/data63/benoitl/Euclid/full_sky/'+suff2+'/'

f=open(path_tiling+'params_tiling.json','r')
par=json.load( f)
f.close()


kmax=par["N_ra"]*par["N_dec"]


list1=''
N=0
for i in range (1, kmax+1):
    if os.path.isfile(path_mask+'mangle/tile_%d/subsample_pixed_balked.pol'%i):
        N+=1
        list1+=path_tiling+'tiles/tiles_%d_%d.pol'%(kmax, i)+' '

os.system('poly2poly  -m1e-8 %s %s'%(list1,path_mask+'footprint.pol' ) ) 
print "there were %d tiles in this footprint"%N
