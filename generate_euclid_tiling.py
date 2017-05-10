import os
import json
import mangle_utils as mu
from Man2Hpx import m_utils as m_u


## Aims at generating a tiling for Euclid to create the Full Sky Mangle masks.  
## A setup is characterized by ra_min, ra_max, dec_min, dec_max, N_ra, N_dec, + a Mangle pixelisation level
## and a path.




dec_max= 89
dec_min=-89
ra_min=0
ra_max=360

N_ra=50
N_dec=50
pix=10

suff='tiling_V2_pix10'

path='/data63/benoitl/Euclid/full_sky/'+suff+'/'

#### make the directories withing the rootpath
if not(os.path.isdir(path)):
    os.mkdir(path)
    os.mkdir(path+'rects/')
    os.mkdir(path+'pols/')
    os.mkdir(path+'tiles/')


par={}

par["dec_min"]=dec_min
par["dec_max"]=dec_max

par["ra_min"]=ra_min
par["ra_max"]=ra_max

par["N_ra"]=N_ra
par["N_dec"]=N_dec
par["pix"]=pix


f=open(path+'params_tiling.json','w')
json.dump(par, f)
f.close()

k=0
for i in range(N_ra):
    for j in range(N_dec):
        k+=1
        print k
        ra0=ra_min+1.0* (ra_max-ra_min)/ N_ra *(i)
        ra1=ra_min+1.0* (ra_max-ra_min)/ N_ra *(i+1)

        
        dec0=dec_min+1.0* (dec_max-dec_min)/ N_dec *(j)
        dec1=dec_min+1.0* (dec_max-dec_min)/ N_dec *(j+1)

        
        f=open(path+'rects/try_%d.rect'%k, 'w')
        s=mu.c2r(ra0, ra1, dec0, dec1)
        print >>f, "rectangle %d"%k
        print >>f, s[1]
        f.close()
        os.system('poly2poly -ir1  -m1e-8   %s %s'%(path+'rects/try_%d.rect'%k, path+'pols/try_%d.pol'%k))
        
kmax=k
## merge all the pols. 

list1=''
for i  in range(1, (kmax+1)/2):
    list1+=(path+'pols/try_%d.pol'%i)+' '

os.system('poly2poly -m1e-8  %s %s '%(list1, path+'full_sphere_l1_%d.pol'%kmax))

list2=''
for i  in range((kmax+1)/2, (kmax+1)):
    list2+=(path+'pols/try_%d.pol'%i)+' '

os.system('poly2poly -m1e-8 %s  %s %s '%(path+'full_sphere_l1_%d.pol'%kmax, list2, path+'full_sphere_%d.pol'%kmax))


os.system('snap -m1e-8 -a0.0001 -b0.0001 -t0.0001  %s %s '%(path+'full_sphere_%d.pol'%kmax, path+'full_sphere_snapped_%d.pol'%kmax))


#### extracting the tiles
for i in range(1, kmax+1):
    os.system('poly2poly -J%d,%d %s  %s'%(i,i,path+'full_sphere_snapped_%d.pol'%kmax, path+'tiles/tiles_snp_%d_%d.pol'%(kmax, i)))
    os.system('pixelize -m1e-8 -Ps0,%d  %s %s '%(par["pix"],path+'tiles/tiles_snp_%d_%d.pol'%(kmax, i),path+'tiles/tiles_%d_%d.pol'%(kmax, i) ))

#### extracting the tiles
#for i in range(1, kmax):
#    os.system('poly2poly -J%d,%d %s  %s'%(i,i,path+'full_sphere_snapped_pixelized_%d.pol'%kmax, path+'tiles/tiles_%d_%d.pol'%(kmax, i)))









