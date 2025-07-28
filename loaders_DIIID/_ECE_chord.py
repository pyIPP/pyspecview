import MDSplus
from time import time
from IPython import embed
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import cumtrapz

try:
    from .map_equ import equ_map
except:
    from map_equ import equ_map


def ece_los(f_los, z_los, Rgrid, Zgrid, B0, ne, Te, niter=5):

    from scipy.constants import c, e,m_e, epsilon_0


    #relativistic mass increase 
    v=np.sqrt(2*Te*(e/m_e))
    gamma = 1/np.sqrt(1-(v/c)**2)
    m_e = m_e*gamma
 
    #plasma and cyclotron frequencies
    wc = e*B0/m_e
    wp2 = ne*e**2/(m_e*epsilon_0)
    w0 = 2*np.pi*f_los
    nharm = 2

    #index squared for extraordinary wave perpendicular to B
    N = np.sqrt(np.maximum(1-wp2/w0**2*(w0**2-wp2)/(w0**2-wp2-wc**2),0))

     
    #follow approach based on Fermat principle
    #https://physics.stackexchange.com/a/606479


    #assume horisontal LOS with only small refraction -> ds/dt = sqrt((dR/dt)**2+(dZ/dt)**2) ~ dR/dt
    R_los = Rgrid[::-1]  #it starts from LFS 
    Z_los = np.zeros_like(R_los)+z_los
    Nspline = RectBivariateSpline(Zgrid,Rgrid,N, kx=2, ky=2)
    

    #calculate Z_los in 3 iterations
    for i in range(niter):
        dNdZ_los = Nspline.ev(Z_los, R_los, dx=1)
        intdNdz = cumtrapz(dNdZ_los, R_los,initial=0)
        N_los = np.maximum(Nspline.ev(Z_los, R_los), 1e-6) #to avoid issues at n=1 resonance
        Z_los = z_los+cumtrapz(intdNdz/N_los, R_los,initial=0)
    
 
    ##calculare resonance location
    wc_los = RectBivariateSpline(Zgrid,Rgrid, wc*nharm).ev(Z_los, R_los)
    ivalid = np.where(np.cumprod(np.diff(wc_los) > 0))[0][-1]+2
    
    
    R_res = np.interp(w0, wc_los[:ivalid], R_los[:ivalid])
    Z_res = np.interp(w0, wc_los[:ivalid], Z_los[:ivalid])
    
 

    #dNdZ = np.gradient(N, Zgrid, axis=0)
    #plt.pcolor(Rgrid,Zgrid, dNdZ  , vmin =-1, vmax=1, cmap='seismic' )
    #plt.colorbar()
    #plt.contour(Rgrid,Zgrid, w0-wc*nharm , colors='w', levels=[0])
    #plt.axis('equal')
    #plt.ylim(Zgrid[0], Zgrid[-1])
    #plt.plot(R_los, Z_los)
    #plt.plot(R_res, Z_res,'o')
    #plt.show()
  


    return R_los, Z_los, R_res, Z_res




def main():
        
    shot  = 175860
    T = 2.0

    f_los = 120e9 #Hz
    z_los = 0




    #load kinetics profiles
    profs = np.load(f'kin_data_{shot}.npz', allow_pickle=True)

    Te = profs['Te'].item()
    rho = Te['rho']
    it = np.argmin(np.abs(Te['tvec']-T))
    Te = Te['data'][it]

    ne = profs['ne'].item()
    it = np.argmin(np.abs(ne['tvec']-T))
    ne = ne['data'][it]


    #load magnetics field
    mdsserver = 'localhost'
    #mdsserver = 'atlas.gat.com'
    MDSconn = MDSplus.Connection(mdsserver)
    eqm = equ_map(MDSconn)
    eqm.Open(shot, 'EFIT02', exp='D3D')

    Rgrid = np.linspace(eqm.Rmesh[0],eqm.Rmesh[-1], 100) 
    Zgrid = np.linspace(eqm.Zmesh[0],eqm.Zmesh[-1], 150) 

    Br,Bz,Bt = eqm.rz2brzt(Rgrid,Zgrid, t_in=T)
    B0 = np.sqrt(Br**2+Bz**2+Bt**2)[0].T


    #get 2D profiles of Te and ne
    R_in,Z_in = np.meshgrid(Rgrid, Zgrid)
    rhoRZ = eqm.rz2rho(R_in[None],Z_in[None], t_in=T, coord_out='rho_tor')[0]

    Te = np.interp(rhoRZ, rho, Te)
    ne = np.interp(rhoRZ, rho, ne)




    t = time()
    for i in range(100):
        ece_los(f_los, z_los, Rgrid, Zgrid, B0, ne, Te)
    print((time()-t)/100)
    

 
if __name__ == "__main__":
    main()

