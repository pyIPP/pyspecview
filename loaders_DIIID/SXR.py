try:
    from loaders_DIIID.loader import * 
except:
    from loader import * 

import os
from multiprocessing import  Pool
import MDSplus as mds
from time import time as T
#TODO lod only a data slice
#TODO calibration for all diagnostics 
#TODO geometry for all diags
from  IPython import embed

from collections import OrderedDict


def check(shot):
    #fastest check if the shotfile exist
    #BUG 
    status = True

    return status

verbose = False




def xraypaths(shot, toroidal=False):

    if toroidal:

        #------------------------------------------------------------------------------
        # Toroidal SXR for shot >= 91000
        # SX45R1S1-12 (Toroidal X-ray), DIII-D shot > 91000 (since '97)
        # 98/05/11 from Snider. Took from EFIT /link/efit/sxrt97.dat
        # 98/05/19 modified by Snider.
        # 2002/01/21 added new sxr locations for shot >= 108690 per G.Jackson
        #------------------------------------------------------------------------------
        if shot >= 177100:
            #from 2019  
            #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zslit_45= array(( 73.1,)*8+( 78.1,)*4)/100
            rslit_45= array((230.4,)*8+(227.2,)*4)/100
            
            rwall_45 = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            zwall_45 = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_45 = rad2deg(arctan2(zwall_45-zslit_45, rwall_45-rslit_45))
            
            
            rwall_195 = array([1.048,1.163,1.240,1.291,1.406,1.547,1.636,1.713,1.879,1.994,2.135,2.275 ]);
            zwall_195 = array([1.129,1.103,1.141,1.219,1.180,1.064,1.012,1.025,1.012,0.999,0.845,0.509 ]);
            zslit_195 = -0.834
            rslit_195 = 2.378
            
            xangle_195 = rad2deg(arctan2(zwall_195-zslit_195, rwall_195-rslit_195))
            
            
                       
            rxray  = r_[rslit_45, rslit_195*ones(12), rslit_195*ones(12)]*100
            zxray  = r_[zslit_45, zslit_195*ones(12), zslit_195*ones(12)]*100
            xangle = r_[xangle_45,xangle_195,xangle_195]

            

                          #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zxray_new= array(( 73.1,)*8+( 78.1,)*4)
            rxray_new= array((230.4,)*8+(227.2,)*4)
            
            Rfin = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            Zfin = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_new = rad2deg(arctan2(Zfin-zxray_new/100, Rfin-rxray_new/100))
       
            rxray_old = zeros(12) + 233.61
            zxray_old = zeros(12) + 78.62 
            xangle_old = array([ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ])
            
            rxray  = r_[rxray_new, rxray_old, rxray_old]
            zxray  = r_[zxray_new, -zxray_old, -zxray_old]
            xangle = r_[xangle_new,-xangle_old,-xangle_old]
      
                    
            #print 'xxx', zxray
        elif shot >= 168847:
            
            #BUG it is wrong!!!!

            #sxrgeom_plot3: load in geometry for 2017 SXR45U view chords
            zxray_new= array(( 73.1,)*8+( 78.1,)*4)
            rxray_new= array((230.4,)*8+(227.2,)*4)
            
            Rfin = array([2.19 , 1.98 ,  1.53, 1.12  , 1.01 , 1.01 , 1.01 ,1.01, 1.01, 1.01, 1.01, 1.01]);
            Zfin = array([-0.83 ,-1.06 ,-1.24 ,-1.34 ,-1.01 ,-0.65 ,-0.28 ,0.06, 0.20, 0.37, 0.61, 0.81]);

            xangle_new = rad2deg(arctan2(Zfin-zxray_new/100, Rfin-rxray_new/100))
       
            rxray_old = zeros(12) + 233.61
            zxray_old = zeros(12) + 78.62 
            xangle_old = [ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ]
            
            rxray  = r_[rxray_new, rxray_old, rxray_old]
            zxray  = r_[zxray_new, zxray_old, zxray_old]
            xangle = r_[xangle_new,xangle_old,xangle_old]

        elif shot >= 134917:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62 
            xangle = [ 227.94, 231.91, 236.02, 240.03, 244.11, 246.11, 
                248.14, 250.07, 252.04, 253.95, 255.83, 259.47 ]
        elif shot >= 130868:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle=[236.124, 240.224, 241.936, 243.979, 246.014, 248.345, 250.409,
                252.682, 256.106, 259.752, 263.364, 267.078]
        elif shot >= 112250:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle = [227.381,229.443,231.505,233.617,237.791,242.040,
                246.281,250.455,254.484,258.415,262.154,265.926]
        elif shot >= 108690:
            rxray = zeros(12) + 233.61
            zxray = zeros(12) + 78.62
            xangle = [227.401,229.523,231.579,233.633,237.713,241.750,
                245.905,249.886,253.810,257.618,261.196,264.833]
        elif shot >= 91000:
            rxray = zeros(12) + 232.0   
            zxray = zeros(12) + 78.56
        #     xangle = [234.60,236.60,238.60,240.70,244.70,248.70,
        #		252.60,256.40,260.10,263.60,266.40,270.00]
            xangle = [225.60,227.60,229.60,231.70,235.70,239.70,
                243.60,247.40,251.10,254.60,257.40,261.00]
        else:
            raise Exception('Toroidal X-ray is not available for shots earlier than 91000.')
        
        if len(xangle) == 12:
            xangle = r_[(xangle,)*3]
            rxray  = r_[(rxray, )*3]
            zxray  = r_[(zxray, )*3]

    else:
        #------------------------------------------------------------------------------
        # Poloidal SXR
        #------------------------------------------------------------------------------
        if shot < 80744:

            #  Old Xray arrays (0, -2 and 2)
            #  Print, 'Xraypaths: Old arrays for shots < 80744'
            xangle=[    120.,124.,128.,132.,136.,140.,144.,148.,  
                152.,156.,160.,164.,168.,172.,176.,180.,180., 
                184.,188.,192.,196.,200.,204.,208.,212.,216., 
                220.,224.,228.,232.,236.,240., 
                294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, 
                262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, 
                234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, 
                202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5]

            zxray =[ -10.7,-10.7,-10.7,-10.7,-10.7,-10.7, 
                -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, 
                -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, 
                -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, 
                -14.7,-14.7, 
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, 
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, 
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, 
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6]

            rxray =[ 248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, 
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, 
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6,  
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, 
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5] 

        elif shot < 91378:
            # sxr0s1-16, sxr0s17-32, sxr2s1-16, sxr2s17-32, DIII-D shot > 80744
            # 94/07/22 from Snider,modified for droop 1-6-95
            # took from EFIT /link/efit/sxr94.dat
            # print,'Xraypaths: for shots between 80744 and 91378.'


            r_0s1 = zeros(16) + 248.9
            z_0s1 = zeros(16) - 10.7
            theta_0s1 = [123.69,126.68,131.24,135.83,140.40,144.89,147.82,154.86,
                    157.87,160.92,163.98,167.06,170.12,173.15,176.13,179.05]

            r_0s17 = zeros(16) + 248.9
            z_0s17 = zeros(16) - 14.7
            theta_0s17= [180.33,183.22,186.18,189.19,192.24,195.30,198.36,201.40,
                    204.41,210.46,214.38,218.87,223.45,228.06,232.64,235.64]

            r_2s1 = zeros(16) + 197.6
            z_2s1 = zeros(16) + 130.1
            theta_2s1 = [290.17,285.65,281.05,277.98,274.92,271.89,268.91,265.98,
                    258.94,255.94,251.89,249.83,246.76,243.72,239.22,234.85]

            r_2s17 = zeros(16) + 194.5
            z_2s17 = zeros(16) + 132.6
            theta_2s17= [233.85,232.41,230.94,227.96,224.94,221.89,218.81,215.74,
                    212.70,209.68,202.71,199.71,196.72,193.69,190.63,187.56]

            xangle = r_[theta_0s1,theta_0s17,theta_2s1,theta_2s17]
            rxray  = r_[    r_0s1,    r_0s17,    r_2s1,    r_2s17]
            zxray  = r_[    z_0s1,    z_0s17,    z_2s1,    z_2s17]

        elif shot < 127721:
        #   New arrays (P1 and M1) from Robin Snider.   tbt

        #   Ted,  Below are the paths for the new x-ray arrays.  The data theta_side
        #   cooresponds to the m1 array.  i.e. the first element cooresponds to
        #   sx90rm1s1 and the last to sx90rm1s32.  The top array cooresponds to
        #   sx90rp1s1 ect.  The shot of interest is 89364.
        #
        #   The shot of interest should be 9137config8 instead of 89364.  The transition
        #   from old soft x-ray system and point names to the new ones for the 
        #   poloidal array was made on April,'97 after the RDP vent. (Snider 7-31-98)

        #   Print, 'Xraypaths: NEW arrays for shots >= 91378'

            theta_side = [98.8,101.6,104.6,107.5,110.5,113.5,115.      
                    ,116.5,118.,119.5,120.9,122.4,123.8,126.7,129.,131.9, 
                    134.8,139.2,143.7,148.1,152.5,156.8,158.8,163.,168.9 
                    ,174.8,180.7,185.,188.8,194.5,199.,207.9]

            z_side = [-77.65,-77.65,-77.65,-77.65,-77.65,-77.65,   
                    -77.65,-77.65,-77.65,-77.65,-77.65,-77.65,  
                    -77.65,-77.65,-77.65,-77.65,-77.65,-77.65,  
                    -77.65,-77.65,-77.65,-77.65,-81.36,-81.36,-81.36  
                    ,-81.36,-81.36,-81.36,-81.36,-81.36,-81.36,-81.36]

            r_side = [237.1,237.1,237.1,237.1,237.1,237.1,   
                    237.1,237.1,237.1,237.1,237.1,237.1,    
                    237.1,237.1,237.1,237.1,237.1,237.1,    
                    237.1,237.1,237.1,237.1,235.6,235.6,    
                    235.6,235.6,235.6,235.6,235.6,235.6     
                    ,235.6,235.6]

            theta_top = [260.9,258.,255.1,252.1,249.1,246.2,243.2   
                    ,240.2,237.3,234.5,230.7,227.9,223.5,219.1,214.6,  
                    210.2,205.9,203.1,201.4,198.6,195.7,192.8,189.9   
                    ,186.9,184.,179.6,175.3,171.5,168.7,165.8,162.9   
                    ,160.]

            z_top = [78.6,78.6,78.6,78.6,78.6,    
                    78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6,78.6  
                    ,78.6,78.6,78.6,82.3,82.3,82.3,82.3,82.3   
                    ,82.3,82.3,82.3,82.3,   
                    82.3,82.3,82.3,82.3,82.3]

            r_top = [236.9,236.9,236.9,236.9,236.9,236.9,   
                    236.9,236.9,236.9,236.9,236.9,236.9,  
                    236.9,236.9,236.9,236.9,236.9,236.9,  
                    235.3,235.3,235.3,235.3,235.3,235.3   
                    ,235.3,235.3,235.3,235.3,235.3,235.3  
                    ,235.3,235.3]

            xangle = r_[theta_top, theta_side]
            rxray  = r_[   r_top,     r_side]
            zxray  = r_[    z_top,     z_side]
            
            
     
        elif shot!=1e6 :
        # New poloidal SXR geometry for 2007 (E. Hollmann). xangle is defined
        # counterclockwise starting from R axis. Chords are RP1 (1:32), then RM1 (1:32).
        
            xangle = [270.8,267.2,263.3,259.2,254.8,250.3,245.6,240.9, 
                236.1,231.4,226.7,222.3,218.0,214.0,210.2,206.7, 
                217.8,214.3,210.6,206.5,202.3,197.8,193.2,188.4, 
                183.6,178.9,174.2,169.7,165.3,161.2,157.4,153.3, 
                89.2,92.8,96.7,100.8,105.2,109.7,114.4,119.1, 
                123.9,128.6,133.3,137.7,142.0,146.0,149.8,153.3, 
                142.2,145.7,149.4,153.5,157.7,162.2,166.8,171.6, 
                176.3,181.1,185.8,190.3,194.6,198.8,202.6,206.7]

            rxray = [230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 	
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 	
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                230.4,230.4,230.4,230.4,230.4,230.4,230.4,230.4, 
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2, 	
                227.2,227.2,227.2,227.2,227.2,227.2,227.2,227.2]
                
            zxray = [73.1,73.1,73.1,73.1,73.1,73.1,73.1,73.1, 
                73.1,73.1,73.1,73.1,73.1,73.1,73.1,73.1, 
                78.1,78.1,78.1,78.1,78.1,78.1,78.1,78.1, 
                78.1,78.1,78.1,78.1,78.1,78.1,78.1,78.1, 
                -73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1, 
                -73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1,-73.1, 
                -78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1, 
                -78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1,-78.1] 		


            xangle = array(xangle)
            #xangle[16:32],xangle[48:64] = copy(xangle[48:64][::-1]),copy(xangle[16:32][::-1])#NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wrong geometry!!!
            
            
        else:
            
            
            xangle = [266.98,263.79,260.41,256.84,253.10,249.21,245.19,241.08,
                      236.92,232.75,228.63,224.58,220.64,216.86,213.24,209.80,
                      215.55,212.18,208.64,204.95,201.11,197.18,193.16,189.11,
                      185.07,181.07,177.15,173.35,169.69,166.19,162.86,159.72,
                      91.03,94.45,98.08,101.90,105.89,110.03,114.28,118.59,
                      122.92,127.21,131.43,135.52,139.46,143.21,146.77,150.12,
                      146.00,149.36,152.91,156.64,160.52,164.54,168.66,172.84,
                      177.02,181.18,185.27,189.24,193.07,196.74,200.22,203.50]
            
            
            rxray = array([230.4]*16+[227.2]*16+[230.4]*16+[227.2]*16)
            zxray = array([ 73.1]*16+[ 78.1]*16+[-73.1]*16+[-78.1]*16)

  

    rxray = array(rxray)/100.
    zxray = array(zxray)/100.

    angle =  deg2rad(xangle) 
    
    
    rxray2 = rxray  + 3 * cos(angle)
    zxray2 = zxray  + 3 * sin(angle)
 
    return rxray,zxray,rxray2,zxray2,angle






shot_year = (
    (0, 2007),
    (127330, 2007),
    (131060, 2008),
    (135130, 2009),
    (140840, 2010),
    (143710, 2011),
    (148160, 2012),
    (152130, 2013),
    (156200, 2014),
    (160940, 2015),
    (164780, 2016),
    (168440, 2017),
    (174580, 2018))

def get_filter(shot, PA):
    
    #find a year corresponding to the discharge number, BUG use SQL? 
    year = 0
    for shot_, year_ in shot_year:
        if  shot < shot_: break
        year = year_
        
    #print year
    
    #exit()
        
    if PA and year < 2007:
        raise Exception('No SXR filter information!')
        #OMFITx.End()
        
    if PA and year in (2007,2008):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick moly', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'50 um thick moly', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'12 um diamond','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'15 um diamond','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'50 um thick moly', 'filter':'32 um diamond','pinhole diameter':1.95e3}]
    if PA and year in (2009,2010):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'?', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'?', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'12 um diamond','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'18 um diamond','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'10 um Ti','pinhole diameter':1.95e3}]
    if PA and year == 2011:
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'?', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'?', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Mid Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'25 um Be','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'?', 'filter':'125 um Be','pinhole diameter':1.95e3}]
    if PA and year in (2012,2013):
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.1},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 2 mm slit', 'foil':'25 um thick moly', 'filter':'none','pinhole diameter':870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'1.3 mm thick tantalum', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Visible light', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'1 mm fused silica','pinhole diameter':1.95e3},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'125 um Be','pinhole diameter':1.95e3}]
    if PA and year in r_[2014:2019]:
        settings=[{'description':'Closed', 'aperture':'closed', 'foil':'none', 'filter':'none','pinhole diameter':0.0},
                  {'description':'Disruption prad', 'aperture':'200 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':200.0},
                  {'description':'ELM prad', 'aperture':'300 um x 3 mm slit', 'foil':'25 um thick moly', 'filter':'none','pinhole diameter':1070.},#870.0},
                  {'description':'Low Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'1.3 mm thick tantalum', 'filter':'12 um Be','pinhole diameter':1.95e3},
                  {'description':'Intermediate prad', 'aperture':'400 um pinhole', 'foil':'50 um thick SS', 'filter':'none','pinhole diameter':400.0},
                  {'description':'Hi Te SXR', 'aperture':'1 mm x 3 mm slit', 'foil':'25 um thick moly', 'filter':'125 um Be','pinhole diameter':1.95e3}]
        
        
        
    if not PA and shot < 168847:
        settings=[{'description':'140um Be', 'aperture':' '     ,'foil':'none',     'filter':'140um Be','pinhole diameter':3.91e3}]

    elif not PA and shot < 172691:
        settings=[
                {'description':'Closed',   'aperture':'closed','foil':'none',     'filter':'none','pinhole diameter':0.0},
                {'description':'10 um Ti', 'aperture':'1 mm x 1.2 mm slit '     ,'foil':'10 um Ti', 'filter':'none','pinhole diameter':3.91e3},
                {'description':'ND1',      'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'none','pinhole diameter':3.91e3},
                {'description':'ND2',      'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'none','pinhole diameter':3.91e3},
                {'description':'25um Be',  'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'25um Be','pinhole diameter':3.91e3},
                {'description':'125um Be', 'aperture':'1 mm x 1.2 mm slit '     ,'foil':'none',     'filter':'125um Be','pinhole diameter':3.91e3}]

        
    elif not PA:
        settings=[
                {'description':'Closed',          'aperture':'closed','foil':'none',    'filter':'none',     'pinhole diameter':0.0},
                {'description':'dummy for bakes', 'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'100 um SS','pinhole diameter':3.91e3},
                {'description':'default SXR',     'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'125 um Be','pinhole diameter':3.91e3},
                {'description':'low Te SXR',      'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'25um Be',  'pinhole diameter':3.91e3},
                {'description':'disruptions Prad','aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'none',     'pinhole diameter':3.91e3},
                {'description':'ELM Prad',        'aperture':'1 mm x 1.2 mm slit ',     'foil':'none',    'filter':'none',     'pinhole diameter':3.91e3}]

                

    return settings

#-----------------------------
# SXR Calibration to turn Volts into W/cm**2
#-----------------------------
def get_calib(shot,calib_path,cam):
    # Get the filter settings options for our shot
    if verbose: print('Getting SXR Filter Settings')
    
    PA = cam in ('90RM1','90RP1')
    filter_settings = get_filter(shot,PA)
    # Get electronic settings and filter setting

    if PA:
        SXRsettings = calib_path+os.sep+'SXRsettingsPA.dat'  
    else:
        SXRsettings = calib_path+os.sep+'SXRsettings45U.dat'  #BUG!!!
    index = open(SXRsettings,'r')

    shots, Rc, Gain,Filt = [],[],[],[]
    do_read = False
    for line in index:
        if do_read:
            P = cam == '90RP1'
            line = line.split()
            shots.append(int(line[0]))
            
            Rc.append(float(line[1+P].replace('k','e3')))
            Gain.append(float(line[2+PA+P]))
            Filt.append(int(line[3+PA*2+P]))
        if line[:4] == 'shot':
            do_read = True
    index.close()
    # Find settings for our shot
    wh = where(array(shots) >= shot)[0]
    # Go back one for the change shot
    try:
        whs = wh[0]-1
    except Exception:
        print('Calibration shot is last shot!')
        whs = len(shots)-1
        
 
    Rc, Gain, filter = Rc[whs], Gain[whs], Filt[whs]
    

    if (shot < 168847 and cam == '45R1') or cam in ('165R1','195R1'):
        #old TA diagnostics
        Rc, Gain, pinhole,filt   = 500e3, 1,  3.91e3,1
    else:
        # Get pinhole diameter (um)
        pinhole = filter_settings[filter]['pinhole diameter']


    if pinhole == 0:
        print('Pinhole closed '+cam)
        pinhole = 1.95e3

 


    return Rc, Gain, pinhole, filter

def get_calib_fact(shot, geometry_path,  toroidal=False):
    ## --------------------------
    # Get subsystem calibrations
    # --------------------------
    if verbose: print('Getting SXR Calibration Settings')
    calib_dict = {}

    if toroidal:
        cams = '45R1','165R1','195R1'
    else:
        cams = '90RM1','90RP1'


    for cam in cams:
        #BUG from for Tor cameras!!!!!!
        calib_path = '/fusion/projects/diagnostics/sxr/'

        try:
            resistor, gain, pinhole, filter = get_calib(shot,calib_path,cam)
        except:
            #if cam == '195R1': shot = min(shot, 160000) #to prevent loading of calibration for new U45 camera
            resistor, gain, pinhole, filter = get_calib(shot,geometry_path,cam)

        # responsivity of AXUV photodiodes [A/W] from E Hollman
        eff_resp =  0.27
        
        # 50 ohm termination divides by 2
        term = 0.5  
        
        if cam in ('90RM1','90RP1'): #poloidal array
            pol_array = True
            nch = 16
            # Effective pinhole thickness (microns)
            pinhole_thick = 25.#50.0
            # Distance center to pinhole (cm)
            ctr_to_pinhole = 2.95#3.00
            # Element area (m^2)
            el_area = 2 * 5/1e6
            # Aperture area (m^2)
            ap_area = 0.25*pi*(pinhole/1e6)**2
            # Distance from center (cm) in eight steps
            dist_from_ctr = ( arange(-nch/2,nch/2) + 0.5 ) * 0.212 #   % [cm]
            tanpsi = dist_from_ctr / ctr_to_pinhole
            cos4 = (tanpsi**2 + 1.)**(-2.)
            thick_factor = abs(tanpsi) * (-4./pi) * (pinhole_thick/pinhole) + 1.

            #measured for filter 5 (high Te)
            if cam == '90RM1':
                etendue  = array( [2.141,2.338,2.528,2.702,2.882,2.695,2.728,2.985,
                                   2.989,2.971,2.858,2.711,2.538,2.349,2.152,1.900,
                                   2.074,2.331,2.589,2.836,3.056,3.234,3.356,3.41,
                                   3.382,3.282,3.122,2.914,2.675,2.419,2.161,1.91])*1e-8
            if cam == '90RP1':
                etendue  = array( [1.845,2.073,2.302,2.523,2.72,2.671,2.813,3.015,
                                   3.031,2.948,2.811,2.631,2.421,2.195,1.965,1.74,
                                   1.881,2.182,2.425,2.663,2.882,3.069,3.211,3.296,
                                   3.314,3.259,3.142,2.975,2.769,2.539,2.297,2.055])*1e-8
                
            if filter == 3:
                #measured, thick slits for low Te
                if cam == '90RM1':
                    etendue  = array( [1.749,2.18,2.442,2.646,2.754,2.558,2.592,2.783,
                                       2.824,2.833,2.783,2.646,2.504,2.309,2.031,1.765,
                                       2.347,2.646,2.895,3.098,3.244,3.352,3.414,3.447,
                                       3.447,3.352,3.211,3.053,2.866,2.625,2.23,1.661])*1e-8
                if cam == '90RP1':
                    etendue  = array( [1.628,2.068,2.356,2.546,2.648,2.607,2.713,2.847,
                                       2.945,2.859,2.758,2.632,2.498,2.315,2.067,1.645,
                                       1.877,2.076,2.292,2.479,2.658,2.766,2.865,2.931,
                                       2.946,2.905,2.803,2.654,2.492,2.23,1.873,1.279])*1e-8
       
            else:
                #correction for different area of the slit 
                etendue *= pinhole**2/1.95e3**2
       
                
            
        if cam in ('45R1','165R1','195R1'): #toroidal array
            pol_array = False
            nch = 20
            # Effective pinhole thickness (microns)
            pinhole_thick = 0.
            # Distance center to pinhole (cm)
            ctr_to_pinhole = 2.7
            # Element area (m^2)
            el_area = 4.1*0.75/1e6
            # Aperture area (m^2)
            ap_area = 0.25*pi*(pinhole/1e6)**2# 12*1/1e6
            #pinhole = ap_area*
            # Distance from center (cm) in eight steps
            dist_from_ctr = ( arange(-nch//2,nch//2) + 0.5 ) * 0.095 #   % [cm]
            tanpsi = dist_from_ctr / ctr_to_pinhole
            cos4 = (tanpsi**2 + 1.)**(-2.)
            thick_factor = abs(tanpsi) * (-4./pi) * (pinhole_thick/pinhole) + 1.

            # overall factor to account for optics...units: m-2 sr-1  it is 1/etendue
            etendue = cos4*thick_factor*el_area*ap_area/ctr_to_pinhole**2*1e4
          
        # Compute calibration
               
        # From V -> W/m^2
        
        #0.35 A/W    U = IR 
        #print 0.35*5e-4*(resistor*gain)
        

        calib = 4*pi/term/(resistor*gain)/eff_resp/etendue
 
        #### Done computing calibration ####
        if pol_array:
            calib_dict[cam+'a'] = calib[:16]
            calib_dict[cam+'b'] = calib[16:]
        else:
            active = [18,17,15,13,11,10,9,8,7,6,5,3]  #from Eric Hollmann

            calib_dict[cam] = calib[active]

        
    return  calib_dict



 
def mds_load(tmp):
    mds_server,  TDI = tmp
    MDSconn = mds.Connection(mds_server )
    data = []
    for tdi in TDI:
        try:
            data.append(MDSconn.get(tdi).data())
        except:
            data.append([])
    return data


def mds_par_load(MDSconn,   TDI,  numTasks=8):
    mds_server = MDSconn.hostspec

    if len(TDI)  == 1:
        return [MDSconn.get(TDI[0])]
        

    #load a junks of a single vector
    TDI = array_split(TDI, min(numTasks, len(TDI)))

    args = [(mds_server,  tdi) for tdi in TDI]

    pool = Pool(len(args))
    out = pool.map(mds_load,args)
    
    
    out = []
    pool = Pool()
    for o in pool.map(mds_load,args):
        out.extend(o)
        
    pool.close()
    pool.join()

    return  out
    



class loader_SXR(loader):
    
    tor_mode_num = True
    pol_mode_num = False  #not implemented yet!! 

    mode_range = (-3,4)
    radial_profile=True
    units = 'W/m$^2$'


    def __init__(self,*args, **kargs):
        super(loader_SXR,self).__init__(*args, **kargs)

   
        self.tree = 'spectroscopy'
        self.phase_diags = '45R1','195R1'
 
        
        self.names = OrderedDict()
        self.names['90RP1'] = self.names['90RM1'] = arange(1,33)
        self.names['45R1' ] = self.names['195R1'] = arange(1,13)

        

        
        if self.shot > 168847: #only one toroidal camera is operational
            self.tor_mode_num = False
        

        self.groups = sorted(list(self.names.keys()))
        path = os.path.dirname(os.path.realpath(__file__))
 
        
        self.calib = get_calib_fact(self.shot,path, True)
        self.calib.update(get_calib_fact(self.shot,path, False))

                
        self.calib['90RM1'] = np.hstack((self.calib.pop('90RM1a'),self.calib.pop('90RM1b')))
        self.calib['90RP1'] = np.hstack((self.calib.pop('90RP1a'),self.calib.pop('90RP1b')))
        #self.calib['90RM1'][2] = 0 #corrupted channel
 
        # Geometry from
        # /usc-data/c/idl/source/efitviewdiagnoses/DIII-D/xraypaths.pro
 
        r_p,z_p,r2_p,z2_p,xangle_p = xraypaths(self.shot, toroidal=False)
        r_t,z_t,r2_t,z2_t,xangle_t = xraypaths(self.shot, toroidal=True)

       
        
        self.Phi = {'165R1':165, '195R1':195, '45R1':45, '90RM1':90, '90RP1':90 }           
        self.R_start = {"90RM1":r_p[32:],"90RP1":r_p[:32],
                       '45R1':r_t[0:12],'165R1':r_t[12:24],'195R1':r_t[24:36]}
        self.R_end =   {"90RM1":r2_p[32:],"90RP1":r2_p[:32],
                       '45R1':r2_t[0:12],'165R1':r2_t[12:24],'195R1':r2_t[24:36]}
        self.z_start = {"90RM1":z_p[32:],"90RP1":z_p[:32],
                       '45R1':z_t[0:12],'165R1':z_t[12:24],'195R1':z_t[24:36]}
        self.z_end =   {"90RM1":z2_p[32:],"90RP1":z2_p[:32],
                       '45R1':z2_t[0:12],'165R1':z2_t[12:24],'195R1':z2_t[24:36]}
        self.theta =   {"90RM1":xangle_p[32:],"90RP1":xangle_p[:32],
                       '45R1':xangle_t[0:12],'165R1':xangle_t[12:24],'195R1':xangle_t[24:36]}
        
 

        self.time_header = {}
        time_header = self.MDSconn.get(f'PTHEAD2("SX90MF01",{self.shot}); __real64').data()[2:]
        self.time_header["90RM1"] = self.time_header["90RP1"] = time_header.reshape(-1,2).T    
        time_header = self.MDSconn.get(f'PTHEAD2("SX195F01",{self.shot}); __real64').data()[2:]
        self.time_header["45R1"]  = self.time_header["195R1"] = time_header.reshape(-1,2).T    

        
        for group, timeh in self.time_header.items():
            if all(timeh==0) and group in self.groups:
                self.groups.remove(group)
        
        if len(self.groups) == 0:
            raise Exception('Fast SXR data are not availible')

        
        self.cache = {cam: {} for cam in self.names.keys()}
        


        
    def get_names(self,group):
        return self.names[group]
                                                    
        
        
    def get_signal(self,groups, names,calib=False,tmin=None,tmax=None):
        
                
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
     
        if size(names) == 1 and not isinstance(names,tuple):
            names = (names,)
       
        if groups == 'all':
            groups = list(self.names.keys())

        if isinstance(groups,str):
            groups = (groups,)
                        
            
            
        idx = []
        TDI = []     
        for group in groups:
            if len(names) == 0:
                channels = self.names[group]
            else:
                channels = np.atleast_1d(np.squeeze(np.int_(names)))
        
            group_ = group.split('R')        
            group_ = group_[0]+group_[1][:-1]
            
            

            indmin = np.where(self.time_header[group][1] > tmin)[0]
            if len(indmin) == 0:
                print(f'tmin {tmin}s is outside of data range {self.time_header[group][1].max()}')
                continue
            else:
                indmin = indmin[0]
            indmax = np.where(self.time_header[group][0] < tmax)[0][-1]+1

            index = np.arange(indmin,indmax)
            #need the first time slice for  background substraction
            if calib: 
                index = unique(r_[0,index])
            
    
            for ch in channels:
                for i in index:
                    name = f'{ch:02d}_{i}'
                    if len(self.cache[group].get(name,[])) <= 1:
                        idx.append([group, ch, i])
                        TDI.append(f'PTDATA2("SX{group_}F{name}",{self.shot},1)')
       


        if len(TDI) > 0:
            out = mds_par_load(self.MDSconn, TDI)
            for [g,ch,i], o in zip(idx, out):
                self.cache[g][f'{ch:02d}_{i}'] = o
  
        #collect all data
        outputs = []


        for g in groups:

            if len(names) == 0:
                channels = self.names[g]
            else:
                channels = np.atleast_1d(np.squeeze(np.int_(names)))
            
            indmin = np.where(self.time_header[g][1] > tmin)[0][0]
            indmax = np.where(self.time_header[g][0] < tmax)[0][-1]+1

            

                
            #collect all data
            data = []
            for ch in channels:
                sig = [self.cache[g][f'{ch:02d}_{i}'] for i in np.arange(indmin,indmax)]
                if  any([len(s) <= 1  for s in sig]):
                    print('data for SXR {g}{ch} are not availible')
                    sig = None
                elif len(sig) == 1:
                    sig = sig[0]
                else:
                    sig = np.hstack(sig)
                data.append(sig)
   
            #prepare time vectors
            nt = max([len(o) for o in data])
            tvec = np.linspace(self.time_header[g][0,indmin], 
                              self.time_header[g][1,indmax-1], nt)
            imin,imax = tvec.searchsorted([tmin,tmax])
            ind = slice(imin,imax+1)
            
            #calibrate data
            if calib:
                data_offset = self.get_signal(g, channels ,calib=False,tmin=-infty,tmax=0)
                if len(names) == 1:
                    data_offset = [data_offset]

                for i,n in enumerate(channels):
                    if len(data[i]):
                        data[i] = data[i]-data_offset[i][1].mean()
                        #if g in self.calib:
                        data[i] *= self.calib[g][int(n)-1]
                
        
            output = []    
            for sxr in data:
                if len(sxr):
                    sxr = sxr[ind]
                else:
                    sxr = np.zeros(imax-imin, dtype='single')
                output.append([tvec[ind], sxr]) 

            outputs += self.hardcoded_corrections(output, g, channels, True)

        if len(output) == 1:
            return outputs[0]
        else:
            return outputs
            
    
                 
   
    def hardcoded_corrections(self,output, camera, channels, fast_data):
        channels = list(channels)
        
        if self.shot < 166887 and camera == '90RP1' and 16 in channels:
            output[channels.index(16)] *= 1.1
            
            
        if self.shot < 166887 and not fast_data and camera == '90RM1' and 10 in channels and 11 in channels:
            ch1 = channels.index(10) 
            ch2 = channels.index(11 ) 
            output[ch1], output[ch2] = output[ch2], output[ch1] 
      
            
        #signal is inversed? 
        if camera == '90RP1' and 22 in channels:
            output[channels.index(22)][1] *= sign(output[channels.index(22)][1])
            
            
        
        if camera == '45R1' and 1 in channels and self.shot  < 180000:
            #BUG od 166566 je to uz vetsi korekce  167451 to zase sedi? 
            output[channels.index(1)][1]*= 1.17
             
            
            
            
        if camera == '195R1':
            if self.shot < 168847 and 1 in channels:
                output[channels.index(1)][1] *= 2*0.96  #od 169350 je to horsi, zeby uz od 168847, od 169593zase prestrelene
            elif self.shot < 175669 and 1 in channels:# 169593:
                output[channels.index(1)][1] *= 2*0.96/1.53
                
            elif self.shot <  181100 and 1 in channels:
                output[channels.index(1)][1] *= 2*0.96/1.53/1.1
            elif 1 in channels:
                output[channels.index(1)][1] *= 2*0.96/1.53/1.1*1.1
            if self.shot >  181100 and 2 in channels:
                output[channels.index(2)][1] *=   1.1
            if fast_data and 11 in channels:
                output[channels.index(11)][1] *=   -1#BUG? in shot 172221

        
        return output
    
        
     
        
     
        
      

    def get_signal_phase(self,name,calib=False,tmin=None,tmax=None):
         
        data = [self.get_signal(d, name,calib=calib,tmin=tmin,tmax=tmax) for d in self.phase_diags]
        tvec = data[0][0]

        return tvec, vstack([d[1] for d in data]).T
    
    def get_names_phase(self):
        names = self.get_names(self.phase_diags[0])
    
        return  names
        
    def get_phi_tor(self,name=None):
        return  deg2rad([self.Phi[d] for d in self.phase_diags])
            
            
    def get_rho(self,group,names,time,dR=0,dZ=0):

        R_start = array([self.R_start[group][int(n)-1] for n in names])
        z_start = array([self.z_start[group][int(n)-1] for n in names])
        R_end = array([self.R_end[group][int(n)-1] for n in names])
        z_end = array([self.z_end[group][int(n)-1] for n in names])
        Phi = array([self.Phi[group] for n in names])
        rho_tg,theta_tg,R,Z = super(loader_SXR,self).get_rho(time,R_start,
                                    z_start,Phi,R_end,z_end,Phi,dR=dR, dZ=dZ)
        
        
        if group in ('45R1','165R1','195R1'):
            rho_tg*= -1
            
            
        return rho_tg, theta_tg,R,Z
    
    
    
    def signal_info(self,group,name,time):

        rho_tg = self.get_rho(group,[ name,],time)[0]
        phi = self.Phi[group]

        info = 'ch:'+str(name)+'  rho_tg: %.2f,  Phi: %d deg'%(rho_tg, phi)

        return info
    
    
 

from matplotlib.pylab import *
def main():
    
    mds_server = "localhost"
    mds_server = "atlas.gat.com"

    import MDSplus as mds
    MDSconn = mds.Connection(mds_server )
    
    from map_equ import equ_map
    eqm = equ_map(MDSconn,debug=False)
    eqm.Open(175860,diag='EFIT01' )
    sxr = loader_SXR(175860,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
    #from ECE import loader_ECE
    #ece = loader_ECE(175900,exp='DIII-D',eqm=eqm,rho_lbl='rho_pol',MDSconn=MDSconn)
  
    t = T()
    out = sxr.get_signal( 'all',[], tmin=4.8, tmax = 4.9,calib=True)
    print(T()-t)
    embed()
    
    exit()
    names = sxr.get_names('90RP1')		
    rho_tg, theta_tg,R,Z = sxr.get_rho('90RM1', names,3)
    #plot(rho_tg);show()

    
    #embed()
    #fit = load('fit_t.npz')
    #plot(tvec, data)
    #plot(fit['tvec'], fit['retro_t'][:,10])
    #show()
    
    tomo_local_path = '~/tomography/'
    loc_dir = os.path.dirname(os.path.abspath(__file__))
    tomo_code_path = tomo_local_path


    tomo_code_path  = os.path.expanduser(os.path.expandvars(tomo_code_path))
    tomo_local_path = os.path.expanduser(os.path.expandvars(tomo_local_path))

    sys.path.append(tomo_code_path)

    #load some modules from pytomo 
    from mat_deriv_B import mat_deriv_B
    from shared_modules import  read_config
    from graph_derivation import graph_derivation
    from geom_mat_setting import geom_mat_setting,prune_geo_matrix
    from annulus import get_bd_mat, get_rho_field_mat
    import config

    shot  = 175900
    input_parameters = read_config(tomo_code_path+"tomography.cfg")
    input_parameters['shot'] = shot
    input_parameters['local_path'] = tomo_local_path
    input_parameters['program_path'] = tomo_code_path
    input_parameters['nx'] = 80
    input_parameters['ny'] = 120

    if not hasattr(config, 'wrong_dets_pref'):
        config.wrong_dets_pref = input_parameters['wrong_dets']
    

    ##if tok_lbl == 'DIIID':
    import geometry.DIIID as Tok
    diag = 'SXR fast'
    diag_path = tomo_local_path+ 'geometry/DIIID/SXR'
    tok = Tok(diag, input_parameters, load_data_only=True, only_prepare=True)#BUG 

    tvec2, data2, data_err = tok.get_data(tmin=3,tmax=3.1)
    
    
    embed()

    plot(tvec2, data2[:,0])
    plot(tvec, data)
    show()
    
    
    calib = get_calib_fact(shot, diag_path,False )
    
    

    
    


    exit()
    #data2 = ece.get_signal("",13, tmin=3.04, tmax = 3.11)
    tvec, sig = data    
    N = 1
    plot(tvec[:(len(sig)//N)*N].reshape( -1,N).mean(1), 
          sig[:(len(sig)//N)*N].reshape(-1,N).mean(-1))
    ##tvec, sig = data2
    #plot(tvec,(sig-mean(sig))/std(sig),'--')
    #plot(tvec,(sig-mean(sig))/std(sig))
    xlim(tvec[0],tvec[-1])
    #xlim( 3.0925,3.092637,)
    show()
    
    exit()
    
      
    
    import IPython
    IPython.embed()
    
    data = array([d for t,d in data1])
    tvec = data1[0][0]
    
        



    polarray()
    
    
    filtered_data = copy(data)
    filtered_data[:12]   -= outer( u1[:,0],noise1)
    filtered_data[12:24] -= outer( u2[:,0],noise2)
    filtered_data[24:]   -= outer( u3[:,0],noise3)
    
    
    plot(tvec.reshape( -1,100).mean(1), data[:12].reshape(12, -1,100).mean(2).T-mean(data[:12,tvec<.0],1))
    
    
    plot(data[3].reshape( -1,100).mean(1).T-mean(data[3,tvec<.1]))
    show()
    
    plot(filtered_data[4].reshape( -1,100).mean(1).T)
    plot(data[4].reshape( -1,100).mean(1).T)
    plot(data[4].reshape( -1,100).mean(1).T- data[-1].reshape( -1,100).mean(1).T)

    show()
    
    plot(s1);plot(s2);plot(s3);show()
    plot(u1[:,0]);plot(u2[:,0]);plot(u3[:,0]);show()

    
    
    imshow(filtered_data[:,:1000], interpolation='nearest', aspect='auto');show()

    imshow(filtered_data.reshape(36, -1,1000).mean(2).T, interpolation='nearest', aspect='auto');show()

    
    plot(data.reshape(36, -1,1000).mean(2).T)
    plot(data[3])
    
    #data1 = sxr.get_signal( '90RP1',range(1,33),tmin=-infty, tmax=infty,calib=False)
    #data2 = sxr.get_signal('90RM1',range(1,33),tmin=-infty, tmax=infty,calib=False)


    data = array([d for t,d in data1]+[d for t,d in data2])#+[d for t,d in data3])
    tvec = data1[0][0]
    
        
    import IPython
    IPython.embed()
    imshow(data[32:].reshape(32, -1,1024).mean(2).T, interpolation='nearest', aspect='auto');show()

    plot(data[:32].reshape(32, -1,256).mean(2).T)

    sig2 = data[32:].reshape(32, -1,4096).mean(2)
    sig2[[2,14]] = 0
    
    
    sig1 = data[:32].reshape(32, -1,1024).mean(2)

    
    plot(sig1.T -sig1[:,:10].mean(1));show()
    imshow(sig2, interpolation='nearest', aspect='auto');show()


    #offset = tvec < .1
    #data_ = data[:,offset]-data[:,offset].mean(1)[:,None]
    #u1,s1,v1 = linalg.svd(data_[:12], full_matrices=False)
    #u2,s2,v2 = linalg.svd(data_[12:24], full_matrices=False)
    #u3,s3,v3 = linalg.svd(data_[24:], full_matrices=False)


    #from scipy import signal
    #fnq = len(tvec)/(tvec[-1]-tvec[0])/2
    #fmax = 50
    #b, a = signal.butter(4, fmax/fnq, 'low')
    #noise1 = inner(u1[:,0], data[:12].T)
    #noise1 -= signal.filtfilt(b,a,noise1)
    #noise2 = inner(u2[:,0], data[12:24].T)
    #noise2 -= signal.filtfilt(b,a,noise2)
    #noise3 = inner(u3[:,0], data[24:].T)
    #noise3 -= signal.filtfilt(b,a,noise3)
    
    
    
    #filtered_data = copy(data)
    #filtered_data[:12]   -= outer( u1[:,0],noise1)
    #filtered_data[12:24] -= outer( u2[:,0],noise2)
    #filtered_data[24:]   -= outer( u3[:,0],noise3)

    #fmax = 3000
    #b, a = signal.butter(6, fmax/fnq, 'low')
    #filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    
    ##fmax = 2000
    ##b, a = signal.butter(6, fmax/fnq, 'low')
    ##filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    #i1 = tvec<.1
    #i2 = tvec> tvec[-1]-.5
    #b1,b2 = filtered_data[:,i1].mean(1), filtered_data[:,i2].mean(1)
    #a1,a2 = tvec[i1].mean(), tvec[i2].mean()

    #filtered_data -= ((b2-b1)/(a2-a1)*(tvec[:,None]-a1)+b1).T
    #offset_err = sqrt(mean(filtered_data[:,tvec<.3]**2,1))
    #error =  offset_err[:,None]+filtered_data*0.05
    #cov_mat = corrcoef(filtered_data[:,tvec<.3])

    
    #f,ax=subplots(2,1,sharex=True, sharey=True)
    #ax[0].plot(tvec, data[:12].T)    
    ##ax[1].plot(tvec, filtered_data.T)
    #ax[1].plot(tvec, filtered_data[:12].T)

    #ax[1].set_xlabel('time [s]')
    #ax[1].set_ylabel('filtered SXR')
    #ax[0].set_ylabel('raw SXR')

    #show()
    
    
    #[errorbar(range(36), d,e) for d,e in zip(filtered_data[:,::10000].T,error[:,::10000].T)];show()

    

    #imshow(filtered_data,interpolation='nearest',aspect='auto',vmax=-0.1,vmin=0.1);colorbar();show()

    
    import IPython
    IPython.embed()
    
    #exit()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #plot(range(28),data[:28,35000:35100])
    #plot(range(28,28*2),data[28:,35000:35100])
    #show()
    
    

    data[30] = 0
    ioff = tvec.searchsorted(tvec[-1]-.1)

    data-= data[:,ioff:].mean(1)[:,None]
    plot(data.mean(1));show()

    
    offset = tvec < .1
    data_ = data[:,offset]-data[:,offset].mean(1)[:,None]
    u1,s1,v1 = linalg.svd(data_[:28], full_matrices=False)
    u2,s2,v2 = linalg.svd(data_[28:], full_matrices=False)

    #plot(range(28),u1[:,0] )
    #plot(range(28,28*2), u2[:,0])
    from scipy import signal
    fnq = len(tvec)/(tvec[-1]-tvec[0])/2
    fmax = 50
    b, a = signal.butter(4, fmax/fnq, 'low')
    noise1 = inner(u1[:,0], data[:28].T)
    noise1 -= signal.filtfilt(b,a,noise1)
    noise2 = inner(u2[:,0], data[28:].T)
    noise2 -= signal.filtfilt(b,a,noise2)
    filtered_data = copy(data)
    filtered_data[:28]  -= outer( u1[:,0],noise1)
    filtered_data[28:]  -= outer( u2[:,0],noise2)
    
    
    fmax = 2000
    b, a = signal.butter(6, fmax/fnq, 'low')
    filtered_data = signal.filtfilt(b,a,filtered_data,axis=1)
    
    
    i1 = tvec<.1
    i2 = tvec> tvec[-1]-.5
    b1,b2 = filtered_data[:,i1].mean(1), filtered_data[:,i2].mean(1)
    a1,a2 = tvec[i1].mean(), tvec[i2].mean()

    filtered_data -= ((b2-b1)/(a2-a1)*(tvec[:,None]-a1)+b1).T
    offset_err = sqrt(mean(filtered_data[:,tvec<.3]**2,1))
    error =  offset_err[:,None]+filtered_data*0.05
    cov_mat = corrcoef(filtered_data[:,tvec<.3])

    #plot(v1[0])
    #plot(v2[0])

    #print data.shape, tvec.shape
    import IPython
    IPython.embed()
    #plot(abs(fft.rfft(data_[0])))

    #NOTE odecist pozadi pred aplikaci IIR  filtru!
    
    f,ax=subplots(2,1,sharex=True, sharey=True)
    ax[0].plot(tvec, data.T)    
    #ax[1].plot(tvec, filtered_data.T)
    ax[1].plot(tvec, filtered_data.T)

    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('filtered SXR')
    ax[0].set_ylabel('raw SXR')

    show()
    
    
    
    offset = filtered_data[:,80000:].mean(1)[:,None]

    contourf(data,20)
    
    #filtered_data-= filtered_data[:,80000:].mean(1)[:,None]
    
    imshow(filtered_data,interpolation='nearest',aspect='auto',vmin=-1000, vmax=1000, extent=(tvec[0], tvec[-1], 0,1));colorbar();show()
    
    
    [errorbar(list(range(64)), d,e) for d,e in zip(filtered_data[:,::10000].T,error[:,::10000].T)];show()

    imshow(filtered_data, interpolation='nearest',aspect='auto');colorbar();show()
    plot((filtered_data-offset)[:,40000]);show()

    import matplotlib.pylab as plt
    plt.plot(tvec, sig)
    plt.show()
    
if __name__ == "__main__":
    main()
    
