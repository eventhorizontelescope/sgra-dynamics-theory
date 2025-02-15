'''  Standard hdf5 loading code. Takes hdf5 files of standard ipole format, and produces 
(1) Iall. An array of shape [x, y, t] of all stokes I images for a model. 
(2) Various useful quantities, e.g. the field of view, flux scale, pixel size, time.
(3) a cylinder plot, q, which shows brightness temperature along the ring as a function of PA and time. 

You may need to change a few quantities for your needsk including
-the for loop values, if youre looking for specific models or another library (v3 vs v5)
-the list_directory, where the files are stored
-A few other changes, e.g. uncommenting out the cylinder plot section, or setting a save_directory or plot_directory in your home folder.

Credit: Based on cylinder.py. Nicholas Conroy, Charles Gammie, Michi Baubock. 

'''

import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
# import scipy.ndimage as ndimage
# import numpy.fft as fft

def myinterp(dat, x, y): ## interpolation
    i = int(round(x-0.5))
    j = int(round(y-0.5))
    ### bilinear:
    di = x-i
    dj = y-j
    z00 = dat[i,j]
    z01 = dat[i,j+1]
    z10 = dat[i+1,j]
    z11 = dat[i+1,j+1]
    idat = z00*(1.-di)*(1.-dj) + z01*(1.-di)*dj + z10*di*(1.-dj) + z11*di*dj
    return idat

###
### Run over all parameters
###

## for v5:
# for flux in [0, 1]: ### 0 is sane, 1 is mad
#     for spin in [-0.97, -0.9375, -0.85, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 0.85, 0.9375, 0.97]: ## double check if it's labeled 0.9375 or 0.94
#         for Rhigh in [40, 160]:
#             for inclination in [30, 50, 90]:

## for v3:
# for window in [3, 4, 5]:
for window in [3]: ### we have windows 3, 4, 5 in v3 library. For whatever reason, I find it works better to run windows separately. 
    for flux in [0, 1]: ### 0 is sane, 1 is mad
        for spin in [-0.94, -0.5, 0,  0.5,  0.94 ]:
            for Rhigh in [1, 10, 40, 160]:
                for inclination in [10, 30, 50, 70, 90, 110, 130, 150, 170]:

                    ###
                    ### get frames and basic parameters
                    ###

                    ### define list directory for v5
                    # list_directory = '/bd6/eht/incite_images/base/0/'
                    # if flux == 1:
                    #     list_directory += 'Ma'
                    # else:
                    #     list_directory += 'Sa'
                    # list_directory += '{0:g}'.format(spin)

                    ### define list directory for v3
                    list_directory = '/bd6/eht/Illinois_SgrA_v3check/230GHz/'
                    if flux == 1:
                        list_directory += 'Ma'
                    else:
                        list_directory += 'Sa'
                    if spin > 0:
                        list_directory += '+'
                    list_directory += '{0:g}_w{1:n}'.format(spin, window)

                    ### load the files
                    files = [f for f in os.listdir(list_directory) if ('h5' in f and f'Rh{Rhigh}' in f and f'i{inclination}' in f)  ] 
                    # files = [f for f in os.listdir(list_directory) if ('h5' in f)] 

                    nftot = len(files)               # total file number
                    sfiles = sorted(files)[:nftot]   # sort because OS returns seemingly random order

                    ### in case we're running the script in a different folder from the files, we must append list_directory pathname to each file name in sfiles
                    slash='/'
                    full_file_names = []
                    for i in np.arange(0, len(sfiles), 1):
                        file_names = list_directory+slash+sfiles[i]
                        full_file_names.append(file_names)
                    sfiles = full_file_names

                    ### set some parameters for the run, based on first file in sequence
                    hfp = h5py.File(sfiles[0],'r')
                    FOV_M = hfp['header']['camera']['dx'][()] # image scale in M
                    FOV_uas =hfp['header']['camera']['fovx_dsource'][()] # image scale in uas
                    if  hfp['header']['camera']['fovx_dsource'][()] != hfp['header']['camera']['fovy_dsource'][()]:
                        print("Error: need a square image. This image appears to be {0:n} x {1:n} pixels".format(hfp['header']['camera']['fovx_dsource'][()], hfp['header']['camera']['fovy_dsource'][()]))
                        system.exit(0)
                    scale = hfp['header']['scale'][()]        # converts pixel to flux density in Jy
                    imagep = np.copy(hfp['pol']).transpose((1,0,2))
                    N = imagep.shape[0]                       # x-size of array in pixels
                    M = imagep.shape[1]                       # y-size of array in pixels
                    dx = FOV_uas/N                            # x-size of pixel in uas
                    dy = FOV_uas/M                            # y-size of pixel in uas
                    half_FOV = FOV_uas/2.

                    # spin = hfp['fluid_header']['a'][()]
                    # inclination = hfp['header']['camera']['thetacam'][()]

                    ### find dt by reading in header/t [which is in units of M] 
                    filefirst = sfiles[0]
                    filelast = sfiles[-1]
                    hf1 = h5py.File(filefirst, 'r')
                    hf2 = h5py.File(filelast, 'r')

                    dt = ( hf2[('header/t')][()] - hf1[('header/t')][()] )/nftot ## framerate in GMc^-3

                    t0 = hf1[('header/t')][()]  # start time
                    Tmax = nftot*dt             # duration
                    tfinal = t0 + Tmax      # end time
                    
                    hf1.close()
                    hf2.close()

                    ###
                    ### now create array 
                    ###

                    ### Allocate space for 2-D Stokes I image for each file
                    Iall = np.zeros((N, M, nftot))

                    ### Loop over all the files to define Iall, which contains the Stokes I data from every frame
                    n = 0
                    print("Reading files from " + str(list_directory) + "...")
                    for fname in tqdm(sfiles):
                        # read in file
                        hfp = h5py.File(fname,'r')
                        imagep = np.copy(hfp['pol']).transpose((1,0,2))
                        nstokes = 0 # 0–4 correspond to Stokes I, Q, U, V, optical depth tau. We want Intensity I.
                        Iall[:,:,n] = imagep[:,:,nstokes] 
                        hfp.close()
                        n += 1

                    ### Smooth. May take a while to smooth if there's a lot of frames like in the v5 library
                    # fwhm = 20.0 # uas. Can change 
                    # sig = fwhm/(dx *(2*np.sqrt(2*np.log(2)))) ## We have 20 uas FWHM resolution. dx = uas/pixels. so 20/dx is FWHM in pixel units.
                    # sIall = ndimage.gaussian_filter(Iall, sigma=(sig,sig,0))

                    ### Get brightness along critical curve
                    # x0 = N/2 + 2*spin*np.sin(inclination * np.pi/180) *(N/FOV_M) # the final factor is a unit conversion. We want in pixel units
                    # y0= M/2 
                    # r_pix = 1*np.sqrt(27)*(N/FOV_M) 
                    # icirc = y0 + r_pix*np.sin(pa)   # for sampling along critical curve
                    # jcirc = x0 + r_pix*np.cos(pa)
                    # x = (jcirc - N/2.)*dx # for plotting
                    # y = (icirc - M/2.)*dy

                    ### Get unsmoothed cylinder plot
                    # q = np.zeros((nftot,ntheta))
                    # for k in range(nftot):
                    #     for l in range(ntheta):
                    #         q[k,l] = myinterp(Iall[:,:,k], icirc[l], jcirc[l]) # interpolate to get brigtness temperature along the ring
                        
                    ### Get smoothed cylinder plot
                    # qs = np.zeros((nftot,ntheta))
                    # for k in range(nftot):
                    #     for l in range(ntheta):
                    #         qs[k,l] = myinterp(sIall[:,:,k], icirc[l], jcirc[l]) # interpolate to get brigtness temperature along the ring

                    ### Get mean im and plot
                    # mean_im = np.mean(Iall, axis = 2) # May take a while to calculate if there's a lot of frames like in the v5 library

                    # fig = plt.figure(figsize = (5,5))
                    # ax = plt.subplot(111)
                    # im = ax.imshow(mean_im.T, cmap='afmhot', aspect = 'auto', origin = 'lower', extent=[-FOV_uas/2, FOV_uas/2, -FOV_uas/2, FOV_uas/2]) # transpose and set origin to lower when plotting with imshow
                    # plt.plot(x, y, 'b--', alpha=0.6)
                    # plt.axis('scaled')
                    # plt.xlim(-FOV_uas/2, FOV_uas/2)
                    # plt.ylim(-FOV_uas/2, FOV_uas/2)
                    # plt.xlabel(r'$x\, [\mu{\rm as}]$')
                    # plt.ylabel(r'$y\, [\mu{\rm as}]$')
                    # plt.show()
                    # plt.close(fig)