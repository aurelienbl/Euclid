import numpy as nm
import os
import matplotlib.pyplot as plt
import sys





def radius_J(mag_tab):
    ind4_5=nm.where( (mag_tab>=4) * (mag_tab<5    ))
    ind5_6=nm.where( (mag_tab>=5) * (mag_tab<6   ))
    ind6_7=nm.where( (mag_tab>=6) * (mag_tab<7    ))
    ind7_8=nm.where( (mag_tab>=7) * (mag_tab<8    ))
    ind8_9=nm.where( (mag_tab>=8) * (mag_tab<9    ))
    ind9_10=nm.where( (mag_tab>=9) * (mag_tab<10   ) )
    ind10_11=nm.where( (mag_tab>=10) * (mag_tab<11  )  )
    ind11_12=nm.where( (mag_tab>=11) * (mag_tab<12 )  )
    ind12_13=nm.where( (mag_tab>=12) * (mag_tab<13  )  )
    ind13_14=nm.where( (mag_tab>=13) * (mag_tab<14  )  )
    ind14_15=nm.where( (mag_tab>=14) * (mag_tab<15 )   )
    ind15_16=nm.where( (mag_tab>=15) * (mag_tab<16 )   )
    ind16_17=nm.where( (mag_tab>=16) * (mag_tab<17)    )
    ind17_18=nm.where( (mag_tab>=17) * (mag_tab<18 )   )
    ind18_19=nm.where( (mag_tab>=18) * (mag_tab<19  )  )
    ind19_20=nm.where( (mag_tab>=19) * (mag_tab<20 )   )


    rad=mag_tab.copy()*0

    rad[ind4_5]=1990
    rad[ind5_6]=1145
    rad[ind6_7]=645
    rad[ind7_8]=330
    rad[ind8_9]=6
    rad[ind9_10]=4
    rad[ind10_11]=3
    rad[ind11_12]=2
    rad[ind12_13]=2
    rad[ind13_14]=1
    rad[ind14_15]=1
    rad[ind15_16]=1
    rad[ind16_17]=0
    rad[ind17_18]=0
    rad[ind18_19]=0
    rad[ind19_20]=0

    return rad





def radius_VIS(mag_tab):
    ind4_5=nm.where( (mag_tab>=4) * (mag_tab<5    ))
    ind5_6=nm.where( (mag_tab>=5) * (mag_tab<6   ))
    ind6_7=nm.where( (mag_tab>=6) * (mag_tab<7    ))
    ind7_8=nm.where( (mag_tab>=7) * (mag_tab<8    ))
    ind8_9=nm.where( (mag_tab>=8) * (mag_tab<9    ))
    ind9_10=nm.where( (mag_tab>=9) * (mag_tab<10   ) )
    ind10_11=nm.where( (mag_tab>=10) * (mag_tab<11  )  )
    ind11_12=nm.where( (mag_tab>=11) * (mag_tab<12 )  )
    ind12_13=nm.where( (mag_tab>=12) * (mag_tab<13  )  )
    ind13_14=nm.where( (mag_tab>=13) * (mag_tab<14  )  )
    ind14_15=nm.where( (mag_tab>=14) * (mag_tab<15 )   )
    ind15_16=nm.where( (mag_tab>=15) * (mag_tab<16 )   )
    ind16_17=nm.where( (mag_tab>=16) * (mag_tab<17)    )
    ind17_18=nm.where( (mag_tab>=17) * (mag_tab<18 )   )
    ind18_19=nm.where( (mag_tab>=18) * (mag_tab<19  )  )
    ind19_20=nm.where( (mag_tab>=19) * (mag_tab<20 )   )


    rad=mag_tab.copy()*0

    rad[ind4_5]=7340
    rad[ind5_6]=4495
    rad[ind6_7]=2725
    rad[ind7_8]=1590
    rad[ind8_9]=805
    rad[ind9_10]=52
    rad[ind10_11]=37
    rad[ind11_12]=23
    rad[ind12_13]=16
    rad[ind13_14]=9
    rad[ind14_15]=5
    rad[ind15_16]=3
    rad[ind16_17]=2
    rad[ind17_18]=1
    rad[ind18_19]=1
    rad[ind19_20]=0

    return rad











def rameam(rr1,rr2):
    if rr1<180:
        rr1+=360
    if rr2<180:
        rr2+=360
    return  nm.mod((rr1+rr2)/2.0  ,360  )

def rameam2(rr1,rr2):
     return  nm.mod((rr1+rr2)/2.0  ,360  )

def rameam3(rr1,rr2):
        if (nm.abs(rr1-rr2))> 180:
            rr11=nm.min([rr1, rr2])
            rr22=nm.max([rr1, rr2])
            rr11+=360
            return nm.mod((rr11+rr22)/2.0, 360)
        else:
            return nm.mod((rr1+rr2)/2.0, 360)

def c2e(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4):
    def rameam(rr1,rr2):
        if rr1<180:
            rr1+=360
        if rr2<180:
            rr2+=360
        return  nm.mod((rr1+rr2)/2.0  ,360  )


    r1=ra1#nm.mod(ra1, 360)
    d1=dec1

    r2= rameam(ra1,ra2)
    d2= 0.5*(  dec1+dec2   )

    r3=ra2#nm.mod(ra2, 360)
    d3=dec2

    r4= rameam(ra2,ra3)
    d4= 0.5*(  dec2+dec3   )

    r5=ra3#nm.mod(ra3, 360)
    d5=dec3

    r6= rameam(ra3,ra4)
    d6= 0.5*(  dec3+dec4   )

    r7=ra4#nm.mod(ra4, 360)
    d7=dec4

    r8= rameam(ra4,ra1)
    d8= 0.5*(  dec4+dec1   )

    s= "%.9f "*16
    a=(r1,d1,r2,d2,r3,d3,r4,d4,r5,d5,r6,d6,r7,d7,r8,d8)

    return a, s%a

def c2e2(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4):
    def rameam(rr1,rr2):
        if (abs(rr1-rr2))> 180:
            rr11=nm.min([rr1, rr2])
            rr22=nm.max([rr1, rr2])
            rr11+=360
            return nm.mod((rr11+rr22)/2.0, 360)
        else:
            return nm.mod((rr1+rr2)/2.0, 360)
      #  if rr1<180:
      #      rr1+=360
      #  if rr2<180:
      #      rr2+=360
   ##     return  nm.mod((rr1+rr2)/2.0  ,360  )

    def decmean(d1, d2):
        return   180./nm.pi *   nm.arcsin (   0.5 *( nm.sin(d1*nm.pi/180.)+nm.sin(d2*nm.pi/180.)    )) 




    r1=ra1#nm.mod(ra1, 360)
    d1=dec1

    r2= rameam(ra1,ra2)
    d2= decmean( dec1,dec2   )

    r3=ra2#nm.mod(ra2, 360)
    d3=dec2

    r4= rameam(ra2,ra3)
    d4= decmean(  dec2,dec3   )

    r5=ra3#nm.mod(ra3, 360)
    d5=dec3

    r6= rameam(ra3,ra4)
    d6=decmean (  dec3,dec4   )

    r7=ra4#nm.mod(ra4, 360)
    d7=dec4

    r8= rameam(ra4,ra1)
    d8= decmean(  dec4,dec1   )

    s= "%.9f "*16
    a=(r1,d1,r2,d2,r3,d3,r4,d4,r5,d5,r6,d6,r7,d7,r8,d8)

    return a, s%a



def c2r(azmin, azmax, elmin, elmax ):
    s= "%f "*4
    a=(azmin, azmax, elmin, elmax )

    return a, s%a







def c2v(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4):
    s= "%f "*8
    a=(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4)

    return a, s%a

def create_weighted_poly_from_edges(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4, id, w, polyname):
    
    
    g=open('jpol', 'w')
    aa=c2e(ra1, dec1, ra2, dec2, ra3, dec3, ra4,dec4)
    print >>g, "edges %d"%id
    print >>g, aa[1]
    g.close()
    os.system('poly2poly -q  -ie2 jpol jpol1')
    
    os.system(' echo %f >> jw'%w)
    os.system(' weight  -q  -zjw jpol1 %s'%polyname)
    os.system(' rm jpol jpol1 jw')
    



