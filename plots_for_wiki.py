import healpy as hp
import matplotlib.pyplot as plt

pathM='/data63/benoitl/Euclid/full_sky/mask_V2/'


coord='EQU'
C='C'



Mm=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_%s_NESTED.fits.gz'%coord, nest=True)
hp.mollview(Mm, nest=True, coord=C, title='Nobs mean o.4096_t.16384')
plt.savefig('./plots_Wiki/full_Mn_%s.png'%coord)
hp.gnomview(Mm, rot=[111, 62.6, 0], nest=True,  coord=C,reso=.05, xsize=2000, title='Nobs mean o.4096_t.16384')
plt.savefig('./plots_Wiki/zoom2_Mn_%s.png'%coord)


M2=hp.read_map(pathM+'/VIS_stars_Nobs_frac2_o.4096_t.16384_%s_NESTED.fits.gz'%coord, nest=True)
hp.mollview(M2, nest=True, coord=C, title='frac2 o.4096_t.16384')
plt.savefig('./plots_Wiki/full_M2_%s.png'%coord)
hp.gnomview(M2, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac2 o.4096_t.16384', coord=C)
plt.savefig('./plots_Wiki/zoom2_M2_%s.png'%coord)


M3=hp.read_map(pathM+'/VIS_stars_Nobs_frac3_o.4096_t.16384_%s_NESTED.fits.gz'%coord, nest=True)
hp.mollview(M3, nest=True, coord=C, title='frac3 o.4096_t.16384')
plt.savefig('./plots_Wiki/full_M3_%s.png'%coord)
hp.gnomview(M3, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac3 o.4096_t.16384', coord=C)
plt.savefig('./plots_Wiki/zoom2_M3_%s.png'%coord)


M4=hp.read_map(pathM+'/VIS_stars_Nobs_frac4_o.4096_t.16384_%s_NESTED.fits.gz'%coord, nest=True)
hp.mollview(M4, nest=True, coord=C, title='frac4 o.4096_t.16384')
plt.savefig('./plots_Wiki/full_M4_%s.png'%coord)
hp.gnomview(M4, rot=[111, 62.6, 0], nest=True, reso=.05, xsize=2000, title='frac4 o.4096_t.16384', coord=C)
plt.savefig('./plots_Wiki/zoom2_M4_%s.png'%coord)






#Mm_EQU=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_EQU_NESTED.fits', nest=True)
#Mm_GAL=hp.read_map(pathM+'/VIS_stars_Nobs_mean_o.4096_t.16384_GAL_NESTED.fits', nest=True)









#hp.mollview(Mm_EQU, nest=True, coord='C', title='Nobs mean o.4096_t.16384 EQU')
#plt.savefig('./plots/full_Mn_EQU.png')
#hp.mollview(Mm_GAL, nest=True, coord='G', title='Nobs mean o.4096_t.16384 GAL')
#plt.savefig('./plots/full_Mn_GAL.png')
