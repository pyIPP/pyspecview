from .loader import * 
import os, logging

logger = logging.getLogger('specview.sece')

'''
Created on Dec 10, 2016
@author: sdenk, todstrci
'''

import numpy as np
import aug_sfutils as sf


def check(shot):
    #fastest check if the shotfile exist

    status = False
    
    shotfiles = 'CTA', 'CTC'
    for s in shotfiles:
        path = shot_path+'%d/XX/%s/%d'
        status |= os.path.isfile(path%(shot//10, s, shot))

    path = shot_path+'%d/L0/%s/%d'
    status |= os.path.isfile(path%(shot//10, 'IEC', shot))
    
    return status

#/afs/ipp-garching.mpg.de/augd/shots/3197/LO/IEC/31978
#/afs/ipp-garching.mpg.de/augd/shots/3197/L0/IEC/


class loader_SECE(loader):

    radial_profile=True

    def __init__(self, *args, **kargs):
        
        super(loader_SECE, self).__init__(*args, **kargs)

        self.Resonances = {}


        self.names = {}
        cta = sf.SFREAD('CTA', self.shot)
        if cta.status:
            names = cta.objects + cta.parsets
            self.names['CTA'] = [n for n in names if n[:2] == 'ch'] [:50]  #last channels are useless     
            try:
                self.Resonances['CTA'] = get_cold_resonances_S_ECE(self.shot, "CTA")
            except Exception as e:
                logger.error(str(e))

        ctc = sf.SFREAD('CTC', self.shot)
        if ctc.status:
            names = ctc.objects + ctc.parsets
            self.names['CTC'] = [n for n in names if n[:2] == 'ch'] [:50]  #last channels are useless      
            try:
                self.Resonances['CTC'] = get_cold_resonances_S_ECE(self.shot, "CTC")
            except Exception as e:
                logger.info(str(e))

        iec = sf.SFREAD('IEC', self.shot)
        if iec.status:
            names = iec.objects + iec.parsets
            self.names['IEC'] = [n for n in names if n[:3] == 'ECE'][:6]  #last channels are useless     
            try:
                self.Resonances['IEC'] = get_cold_resonances_S_ECE(self.shot, "IEC")
            except Exception as e:
                logger.info(str(e))
                
        self.groups = list(self.names.keys())

        self.shotfile = ''
        
        logger.info( 'initialized')


    def get_names(self, group):
        return self.names[group]
    
    def get_signal(self, group, name, calib=False, tmin=None, tmax=None):
        
        if isinstance(name, str):
            names = (name, )
        else:  names = name
        
        
        if tmin is None:    tmin = self.tmin
        if tmax is None:    tmax = self.tmax
        
        if not self.shotfile == group:
            sfo = sf.SFREAD(group, self.shot, experiment=self.exp, edition=self.ed)
        tvec = sfo.gettimebase(names[0])
        noff, nbeg, nend = tvec.searchsorted((0, tmin, tmax))
        self.shotfile =  group

        out = []
        for n in names:
            sig = sfo.getobject(n, cal=calib, nbeg=nbeg, nend=nend)
            if calib:
                sig -= sfo.getobject(n, cal=calib, nend=noff+10).mean()
                if group in ('CTA', 'CTC'):  sig*= -1
                
            out.append([tvec[nbeg:nend+1], sig])
  
        if len(names) == 1: return  out[0]
        
        return out


    def signal_info(self, group, name, time):

        rho, theta, R, z = self.get_rho(group, [name, ], time)

        info = 'ch: '+str(name)+'  R:%.3fm   z:%.3fm  '%(R, z)+self.rho_lbl+': %.3f'%rho
        return info
    

    def get_RZ_theta(self, time, group, names, dR=0, dZ=0, cold=False):

        Brt, Bzt, Btt = sf.Bmesh(self.eqm, t_in=time)
        Br, Bz, Bt = Brt[0], Bzt[0], Btt[0]
        B_spl = RectBivariateSpline(self.eqm.Rmesh, self.eqm.Zmesh, np.sqrt(Br**2 + Bt**2 + Bz**2))

        ch_num = []
        for n in names:
            for  i, a in enumerate(self.names[group]):
                if a == n: ch_num.append(i)
        
        logger.info('get resonance positions')

        if group in self.Resonances:
            R, z = self.Resonances[group](time, self.eqm.Rmesh[0], self.eqm.Rmesh[-1], 
                                               self.eqm.Zmesh[0], self.eqm.Zmesh[-1], B_spl, ch_num)
        else:
            R, z = zeros_like(ch_num), zeros_like(ch_num)

        r0 = np.interp(time, self.eqm.time, self.eqm.Rmag) + dR
        z0 = np.interp(time, self.eqm.time, self.eqm.Zmag) + dZ

        return R, z, np.arctan2(z-z0, R-r0)


    def get_rho(self, group, names, time, dR=0, dZ=0):

        logger.info('SECE: get_rho %.2f %.2f %.2f', time, self.tmax, self.tmin)
        time = max(min(time, self.tmax), self.tmin)

        R, z, theta = self.get_RZ_theta(time, group, names, dR=dR, dZ=dZ)

        rho = sf.rz2rho(self.eqm, R-dR, z-dZ, time, coord_out = self.rho_lbl)[0]
        r0 = np.interp(time, self.eqm.time, self.eqm.Rmag) + dR
        R = np.atleast_1d(R)
        rho[R < r0] *= -1

        return rho, theta, R, z


import ctypes as ct
libECRHso = '/afs/ipp/u/ecrh/lib/libaug_ecrh_setmirrors.so'
from scipy.interpolate import InterpolatedUnivariateSpline , RectBivariateSpline
import numpy as np
import datetime
# NOTE: This routine is not compatible with ECRH 3, yet.
# Contact S. Denk for a new version.

gy_pos_x = np.r_[(2.380, )*2, (2.311, )*2, (2.361, )*4]
gy_pos_y = 0., 0., -0.075, 0.075, -0.115, 0.115, 0.115, -0.115
gy_sect  = np.r_[(7., )*4, (4.5, )*4]
gy_pos_z = np.r_[(0., )*4, (0.32025, )*2, (-0.32025, )*2]
gy_curv_y_105 = np.r_[(0.8793, )*2, (2.9664, )*2, (1.1700, )*4 ]
gy_curv_y_140 = np.r_[(0.8793, )*2, (2.9664, )*2, (0.8540, )*4 ]
gy_curv_z_105 = np.r_[(0.8793, )*2, (2.9664, )*2, (1.1700, )*4 ]
gy_curv_z_140 = np.r_[(0.8793, )*2, (2.9664, )*2, (0.8540, )*4 ]
gy_width_y_105= np.r_[(0.0364, )*2, (0.0329, )*2, (0.0229, )*4 ]
gy_width_y_140= np.r_[(0.0364, )*2, (0.0329, )*2, (0.017, )*4 ]
gy_width_z_105= np.r_[(0.0364, )*2, (0.0329, )*2, (0.0229, )*4 ]
gy_width_z_140= np.r_[(0.0364, )*2, (0.0329, )*2, (0.017, )*4 ]


gy_name = np.array(['ECRH1_1', 'ECRH1_2', 'ECRH1_3', 'ECRH1_4', \
                    'ECRH2_1', 'ECRH2_2', 'ECRH2_3', 'ECRH2_4'])


class libECRH_wrapper:
    def __init__(self, shot):
        shot_folder = '/afs/ipp-garching.mpg.de/augd/shots/'
        path = shot_folder + '%d/' + 'L0/' + 'MAI' + '/%d'
        path = path % (shot // 10, shot)
        date_obj = datetime.datetime.fromtimestamp(os.path.getmtime(path))
        self.date = ct.c_double(date_obj.date().year * 1.e4 + date_obj.date().month * 1.e2 + date_obj.date().day + 0.0)

    def setval2tp(self, gy_obj, gynum, a, beta):
        ct_gy = ct.c_int32(gynum)
        error = ct.c_int32(0)
        error.value = 0
        ct_theta_pol = ct.c_double(a)
        ct_phi_tor = ct.c_double(beta)

        try:
            libECRH = ct.cdll.LoadLibrary(libECRHso)
            libECRH.setval2tp_(ct.byref(error), ct.byref(ct_gy), ct.byref(ct_theta_pol), ct.byref(ct_phi_tor), ct.byref(self.date))

        except:
            raise Exception('ECRH library %s was not found'%libECRHso)
        
        if(error.value != 0):
            logger.error('Encountered error %d', error.value)
            if(np.abs(error.value) == 102 or np.abs(error.value) == 2):
                gy_obj.error = -1
        theta_pol = ct_theta_pol.value
        phi_tor = ct_phi_tor.value  # note - phi tor changes during poloidal sweep
        return theta_pol, phi_tor


class gyrotron:

    def __init__(self, shot, N, ECS, ECN, view=False):  # N is the gyrotron number from 1-8 (ECRH III not included)

        self.x = gy_pos_x[N - 1]
        self.y = gy_pos_y[N - 1]
        self.z = gy_pos_z[N - 1]
        sect = gy_sect[N - 1]
        self.name = gy_name[N - 1]
        self.R = np.hypot(self.x,  self.y )
        self.phi = np.arctan2(self.y, self.x) + sect / 8.0 * np.pi
        self.x = np.cos(self.phi) * self.R
        self.y = np.sin(self.phi) * self.R
        libECRH_obj = libECRH_wrapper(shot)
        self.error = 0
        N_2 = (N-1)%4 + 1

        gyro_name = 'P_sy%d_g%d'%((N-1)//4+1, (N-1)%4+1)
        gyro_pset = ECS(gyro_name)
        self.f = gyro_pset['gyr_freq']
        if(np.abs(self.f - 1.40e11) < 5.e9):
            self.curv_y = gy_curv_y_140[N - 1]
            self.curv_z = gy_curv_z_140[N - 1]
            self.width_y = gy_width_z_140[N - 1]
            self.width_z = gy_width_z_140[N - 1]
        elif(np.abs(self.f - 1.05e11) < 5.e9):
            self.curv_y = gy_curv_y_105[N - 1]
            self.curv_z = gy_curv_z_105[N - 1]
            self.width_y = gy_width_z_105[N - 1]
            self.width_z = gy_width_z_105[N - 1]
        elif(self.f < 5.e9 and not view):
            self.avail = False
            logger.warning('Gyrotron %s not available' %self.name)
            self.error = -1
            return
        elif(not view):
            logger.warning('Found gyrotron with freq=%.2f which is currently not supported', self.f)
            self.error = -1
            return
        if self.curv_y == 0.:
            logger.error('Zero encountered in curvature')
            logger.error('Gyrotron data not properly read')
        self.avail = True
        if(N < 4):
            self.a    = gyro_pset['GPolPos'] * 1000
            self.beta = gyro_pset['GTorPos']
        elif(N < 9):
            SUCOMDAT = ECN('SUCOMDAT')
            if(N == 5):
                beta = SUCOMDAT[26] / 10.0
                a = SUCOMDAT[27] / 10.0
            elif(N == 6):
                beta = SUCOMDAT[71] / 10.0
                a = SUCOMDAT[72] / 10.0
            elif(N == 7):
                beta = SUCOMDAT[116] / 10.0
                a = SUCOMDAT[117] / 10.0
            elif(N == 8):
                beta = SUCOMDAT[161] / 10.0
                a = SUCOMDAT[162] / 10.0
            a_ECS    = gyro_pset['GPolPos'] * 1000
            beta_ECS = gyro_pset['GTorPos']
            if(np.abs(a - a_ECS) > 0.1):
                logger.error('ECS vs SUCOMDAT: different spindle position')
                logger.error('ECS: %.2f, SUCOMDAT: %.2f', a, a_ECS)
                raise ValueError
            if(np.abs(beta - beta_ECS) > 0.1):
                logger.error('ECS vs SUCOMDAT: different toroidal launching angle')
                logger.error('ECS: %.2f, SUCOMDAT: %.2f', a, a_ECS)
                raise ValueError

        if(N < 4):
            gy = 100 + N
        elif(N < 9):
            gy = 200 + N - 4
        else:
            logger.warning('Gy > 8 not supported')
            self.error = -1
        if(N < 5):
            self.time = ECS('T-B')
            self.PW = ECS("PG{0:n}".format(N_2)).astype(np.double)
            self.theta_pol, self.phi_tor = libECRH_obj.setval2tp(gy, self.a, self.beta)
        elif(N < 9):
            time = ECS('T-B')

            self.time = ECN('T-Base')
            self.a_t = ECN("G{0:n}POL".format(N_2)) * 10.0

            a_t_0 = self.a_t[0]
            self.a_t = self.a_t - self.a_t[0] + a_ECS  # Treat last a_t as offset and replace it with a from ECS
# According to M. Schubert the offset might also drift
# -> get the value at the beginning and at the end of the discharge
            self.beta_t = np.zeros(len(self.a_t))
# Changes is beta during discharge not supported
            self.beta_t[:] = beta
# Treat inital error as offset beta_t - beta_t[0] +
            PW = np.double(ECS("PG{0:n}N".format(N_2)))
            pw_spline = InterpolatedUnivariateSpline(time, PW)
            self.PW = pw_spline(self.time)
            if view :
                self.avail = True
            else:
                self.avail = np.any(self.PW > 5.e3)
            if not self.avail:
                logger.error('Gyrotron %s not active', self.name)
                return
            self.theta_pol = np.zeros(len(self.time))
            self.phi_tor = np.zeros(len(self.time))
            for i in range(len(self.time)):
                self.theta_pol[i], self.phi_tor[i] = libECRH_obj.setval2tp(self, gy, self.a_t[i], self.beta_t[i])
        else:
            logger.warning('Gy > 8 not supported')


def get_ECRH_viewing_angles(shot, LOS_no):

    ECS = sf.SFREAD("ECS", int(shot))
    ECN = sf.SFREAD("ECN", int(shot))

    return gyrotron(shot, LOS_no, ECS, ECN, True)


def get_freqs(shot, diag):
    if(diag.name == "IEC"):
        f = []
        dfreq_IEC = 3.0e9
        for i in range(6):
            f.append(132.5e9 + i * dfreq_IEC)
        f = np.array(f)
    elif(diag.name == "CTC"):
        f = [137.0000, 137.6500, 138.0750, 138.3750, 138.5700, \
                 138.6600, 138.7400, 138.8200, 138.9000, 138.9800, \
                 139.0600, 139.1400, 139.2200, 139.3000, 139.3800, \
                 139.4600, 139.5400, 139.6200, 139.7000, 139.7800, \
                 139.8600, 139.9400, 140.0200, 140.1000, 140.1800, \
                 140.2600, 140.3400, 140.4200, 140.5000, 140.5800, \
                 140.6600, 140.7400, 140.8200, 140.9000, 140.9800, \
                 141.0600, 141.1400, 141.2800, 141.5300, 141.8800, 142.3550, 143.0000]
        f = np.array(f) * 1.e9
    elif(diag.name == "CTA"):
        f = np.array([ 135.57, 136.32, 136.82, 137.32, 137.82, 138.12, 138.22, \
            138.32, 138.42, 138.52, 138.62, 138.72, 138.82, 138.92, \
            139.02, 139.12, 139.22, 139.32, 139.42, 139.52, 139.62, \
            139.72, 139.82, 139.92, 140.02, 140.12, 140.22, 140.32, \
            140.42, 140.52, 140.62, 140.72, 140.82, 140.92, 141.02, \
            141.12, 141.22, 141.32, 141.42, 141.52, 141.62, 141.72, \
            141.82, 141.92, 142.02, 142.32, 142.82, 143.32, 143.82, \
            144.57]) * 1.e9

    return f


class Diag:

    def __init__(self, name, exp, diag_str, ed):

        self.name = name
        self.exp = exp
        self.diag = diag_str
        self.ed = ed


class get_cold_resonances_S_ECE:

    def __init__(self, shot, diag_name ):

        self.shot = shot
        self.diag_name = diag_name
        
        if (diag_name == "CTC" or diag_name == "IEC"):
            beamline = 5
        elif diag_name == "CTA":
            beamline = 6
        else:
            logger.error('Unknown Diag name: %s', diag_name)
            raise ValueError
            
        self.gy = get_ECRH_viewing_angles(shot, beamline)


    def __call__(self, time, R_min, R_max, z_min, z_max, B_spline, ch_no, exp="AUGD", diag="None", ed=0):
    # B_spline is expected to be a UnivariateSpline of the total magnetic field
        

        import scipy.constants as cnst

        diagnostic = Diag(self.diag_name, exp, diag, ed)
        f = get_freqs(self.shot, diagnostic) 

        x = self.gy.x, self.gy.y, self.gy.z
        norm = np.hypot(self.gy.R , self.gy.z)
        t1 = np.argmin(np.abs(self.gy.time - time + 0.005))
        t2 = np.argmin(np.abs(self.gy.time - time - 0.005))
        if(t1 == t2):   t2 += 1
        
        phi_tor = -np.mean(self.gy.phi_tor[t1:t2]) / 180.0 * np.pi
        theta_pol = -np.mean(self.gy.theta_pol[t1:t2]) / 180.0 * np.pi
        k = 1.0, np.arctan2(-x[1], -x[0]) + phi_tor, np.pi / 2.e0 + theta_pol
        k_x = np.array([np.cos(k[1]) * np.sin(k[2]), np.sin(k[1]) * np.sin(k[2]), np.cos(k[2])])

        los = [x]
        s = [0.0]    
        los_finished = False
        
        while(True):
            next_step = los[-1] + k_x / 10.0
            R_ray = np.hypot(next_step[0], next_step[1])
            if(R_ray > R_max or R_ray < R_min or next_step[2] > z_max or next_step[2] < z_min):
                break  
            
            los.append(next_step)
            s.append(s[-1] + 0.1)
            
        los = np.array(los).T
        s = np.array(s)
        R_ray = np.hypot(los[0],  los[1])
    
        f_cyc_2 = B_spline(R_ray , los[2], grid=False) * cnst.e / cnst.m_e / cnst.pi
        R_spl = InterpolatedUnivariateSpline(s, R_ray,  k=1)
        z_spl = InterpolatedUnivariateSpline(s, los[2], k=1)
        
        ch_no = np.atleast_1d(ch_no)
        R_rez = []
        Z_rez = []
        for ch in ch_no:
            root = InterpolatedUnivariateSpline(s, f_cyc_2 - f[ch - 1]).roots()
            if len(root) == 0:
                iroot = np.argmin(np.abs(f_cyc_2 - f[ch - 1]))
                R_rez.append( R_ray[iroot])
                Z_rez.append(los[2][iroot])
            else:
                R_rez.append(R_spl(root[0]))
                Z_rez.append(z_spl(root[0]))

        return np.array(R_rez), np.array(Z_rez)
