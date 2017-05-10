import healpy as hp
import matplotlib.pyplot as plt

pathM='/data63/benoitl/Euclid/full_sky/mask_V2/'

Mm=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_ECL_NESTED.fits', nest=True)
M16=hp.read_map(pathM+'/VIS_stars_Nobs_tiny_t.16384_ECL_NESTED.fits', nest=True)
M2=hp.read_map(pathM+'/VIS_stars_Nobs_frac2_o.4096_t.16384_ECL_NESTED.fits', nest=True)
M3=hp.read_map(pathM+'/VIS_stars_Nobs_frac3_o.4096_t.16384_ECL_NESTED.fits', nest=True)
M4=hp.read_map(pathM+'/VIS_stars_Nobs_frac4_o.4096_t.16384_ECL_NESTED.fits', nest=True)


Mm_EQU=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_EQU_NESTED.fits', nest=True)
Mm_GAL=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_GAL_NESTED.fits', nest=True)


hp.mollview(Mm, nest=True, coord='E', title='Nobs mean o.4096_t.16384')
plt.savefig('./plots/full_Mn.png')
hp.mollview(M16, nest=True, coord='E', title='Nobs t.16384' )
plt.savefig('./plots/full_M16.png')
hp.mollview(M2, nest=True, coord='E', title='frac2 o.4096_t.16384')
plt.savefig('./plots/full_M2.png')
hp.mollview(M3, nest=True, coord='E', title='frac3 o.4096_t.16384')
plt.savefig('./plots/full_M3.png')
hp.mollview(M4, nest=True, coord='E', title='frac4 o.4096_t.16384')
plt.savefig('./plots/full_M4.png')

hp.gnomview(Mm, rot=[111, 62.6, 0], nest=True, reso=.1, xsize=2000, title='Nobs mean o.4096_t.16384')
plt.savefig('./plots/zoom1_Mn.png')
hp.gnomview(M16, rot=[111, 62.6, 0], nest=True, reso=.1, xsize=2000, title='Nobs t.16384')
plt.savefig('./plots/zoom1_M16.png')
hp.gnomview(M2, rot=[111, 62.6, 0], nest=True, reso=.1, xsize=2000, title='frac2 o.4096_t.16384')
plt.savefig('./plots/zoom1_M2.png')
hp.gnomview(M3, rot=[111, 62.6, 0], nest=True, reso=.1, xsize=2000, title='frac3 o.4096_t.16384')
plt.savefig('./plots/zoom1_M3.png')
hp.gnomview(M4, rot=[111, 62.6, 0], nest=True, reso=.1, xsize=2000, title='frac4 o.4096_t.16384')
plt.savefig('./plots/zoom1_M4.png')


hp.gnomview(Mm, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='Nobs mean o.4096_t.16384')
plt.savefig('./plots/zoom2_Mn.png')
hp.gnomview(M16, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='Nobs t.16384')
plt.savefig('./plots/zoom2_M16.png')
hp.gnomview(M2, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac2 o.4096_t.16384')
plt.savefig('./plots/zoom2_M2.png')
hp.gnomview(M3, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac3 o.4096_t.16384')
plt.savefig('./plots/zoom2_M3.png')
hp.gnomview(M4, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac4 o.4096_t.16384')
plt.savefig('./plots/zoom2_M4.png')


hp.mollview(Mm_EQU, nest=True, coord='C', title='Nobs mean o.4096_t.16384 EQU')
plt.savefig('./plots/full_Mn_EQU.png')
hp.mollview(Mm_GAL, nest=True, coord='G', title='Nobs mean o.4096_t.16384 GAL')
plt.savefig('./plots/full_Mn_GAL.png')
