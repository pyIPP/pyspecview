# Introduction
The control of the GUI is simple very intuitive but to use the full ability of this tool, 
this short manual was prepared:

## Supported Machines
DIII-D, NSTX, ASDEX Upgrade

## How to start the program
    Execute command 'pyspecview' in the linux clusters toks01-06 or toki01-07. 
    Make sure that the server is not heavily loaded by other users. (command htop)

## Spectrogram panel
    1) set shot, diagnostic, group, and signal
    2) zoom by a mouse wheel, pan by  pressing the mouse middle button
      or by a standard matplotlib zooming and panning toolbar 
    3) set optimal time resolution and contrast by a GAMMA and TAU sliders
    4) Method menu gives you the option to switch between various versions 
       of short time fourier transformation.  The Sparse option is excellent
       for a coherent MHD modes with sufficient signal level.     
    extra)  double clicking will trace a mode, still in the development phase 
    
 <img src="https://user-images.githubusercontent.com/32073690/197647117-ffd1acea-816e-40d7-b762-6c959dfd6f57.png"   height="800"  align="middle">

## Cross-phaseogram panel
    1) set discharge and signal 
    2) select time and frequency like in the Spectrogram panel
    3) by pressing a single color in the colorbar, it will set the other mode numbers to white
    4) Elm coils (ballooning coils), and SXR measure the toroidal mode number 
       Mirnov coils the poloidal mode number. But be very careful in the interpretation
       of the poloidal mode numbers! it never works for 1/1 modes and also other 
       modes do not work very well. Probably the best tool to identify 
       the poloidal mode number is a proper SXR tomography. 
    
    
## Radial profile
    1) In the Spectrogram window select a desired frequency and time range by a right mouse button
    2) for the cross-phase measurement is used the signal with the highest mode amplitude
    but alternatively any other diagnostic can be used (SXR is very good for the core)
    3) when everything is fine, the circles and squares should overlap!  
    4) the black crosses in the cross-phase plot correspond to the expected phase-shift for m=1 mode.
    If the measured phase does not match with expectations, correction in the equilibrium
     position should be introduced by dR and dZ in units of meters. 
    (the difference can also be introduced by the fact that no ray-tracing 
    (and cold resonance for ECE resonance are used)
    
 <img src="https://user-images.githubusercontent.com/32073690/197647365-f83b29cf-2c35-4496-8d58-b27a164d99e3.png"   height="800"  align="middle">


## 2D Te
    1) the same procedure as for Radial profile. The phase and amplitude is then
    mapped on the 2D coordinates in a real space
    2) Poloidal mode number must be set by a user! 
    3) Mouse wheel can be used to move in the time 
    4) Scale slider set a lower limit of the colorbar
    5) Substract option switches between the absolute temperature and the mean subtracted value. 
    NOTE 1: Only LFS measurements are used! 
    
 Example of 2D Te profile of 3/2 NTM
 
 <img src="https://user-images.githubusercontent.com/32073690/198147512-2f9136ec-4172-45e4-95a9-4eb182dd85ed.png"   height="800"  align="middle">

 
## 2D SXR
    The tomographic code will reconstruct stationary emissivity and several harmonics of the selected mode.
    1) Time and frequency are defined by the same procedure as for Radial profile ot 2D Te. 
    2) The user must specify M and N number of the mode (N is necessary only when F camera is used or  the mode
        is compared with ECE at different toroidal locations).
    3) Remove the corrupted channels by "Set data" button. It can be done in the same fashion as in pyTOMO
        - by double click on the channel. Or if the discharge was already    reconstructed by pyTOMO, 
          it will use load setting from there. 
    4) Finally, the quality of the fit should be checked by "Retrofit" button. Reconstruction should match
        amplitude and phase very well. If not, check that the M and N number are right (try to change
        their sign?), or check of there are no remaining corrupted channels or the regularisation is too strong
    5) Checkbox Te contours is used to view conrours of constant Te
    6) Background checkbox allowes tp substract time averaged profile
    7) Regularisation slider - set how much will be imposed smoothes constrain
 
    NOTE: If there are two modes locked at the same frequency, it should be still possible to reconstruct them, but I have not tried, it can be tricky.  
 
 Example of SXR profile from 2/2 mode overlapped by black contours of dTe
 <img src="https://user-images.githubusercontent.com/32073690/198327865-eb159c3d-6c95-4eba-965d-5815f774a612.png"   height="800"  align="middle">

   
 
## General comments:
    Plot can be saved by File/Save plot
    All data, spectrogram, profiles, etc. are stored by File/Save data
    

 
    
