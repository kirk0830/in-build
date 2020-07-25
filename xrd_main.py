import math
import numpy as np
import matplotlib.pyplot as plt
#import time

print('\n\n')
print('!-----------------------------------------------------------------------!')
print('!                                                                       !')
print('!                                                                       !')
print('!                      X-ray diffraction simulation                     !')
print('!                                                                       !')
print('!                                                   ykhuang@dicp.ac.cn  !')
print('!-----------------------------------------------------------------------!')

print('\nDeveloper notes:\n'
    +'1. up to now, only single element material is supported. :(\n'
    +'2. supercell function will come soon, but that will be expensive if used.\n'
    +'3. i will seperate the whole program into several subroutine in near \n'
    +'   future, in order to freely add more functions into it.\n'
)

#-------------------------------INPUT SECTION-------------------------------------
# input section, run type control
pxrd = True
print_debug_mode = True # if activated, only print input parameters, no calculation
matrix_init_debug = False # if activated, only print input parameters and initialized matrix, no calculation

# input section, control physical parameters
wvl_inp = 0.1 # incident ray, [nm]
dist_inp = 1 # distance from sample to screen, [m]
angle_max_inp = 45
angle_f = False
angle_res_func = [0, 0, 0, 0] # angle resoluted function, polynomial formatted?

# input section, accuracy parameters
resol = 2.5
# resolution parameter, axis will be divided into 10^resol pieces. An non-integar is not recommended (although it can work).
# resol = 2, cheap, not accurate
#       = 3, expensive, much more accurate than resol = 2
coord_scale = 2 # 1: no scaling, 2: inflate all coordinates to 10 times, 3: 100 times, ...

denoise = 'cutoff' 
# denoise method, 
# cutmean: only preserve values larger than average, 
# cutoff: only preserve number of values, false: do not discard any data
denoise_cutoff = 0.5
# denoise parameter, meaning differs from methods.
# cutoff: only data points larger than cutoff*max will preserve
# cutmean: only data points larger than (1-cutoff)*mean will preserve

smooth = True
# gaussian smooth on 2d-diffraction pattern
smooth_cutoff = 5 # unit [pixels]
# gaussian smooth accelerate parameter, small value will give rough pattern
# if set as 10, gaussian smooth will calculate only from x[xi-10, xi+10], y[yi-10, yi+10], 
# i.e., a square-shaped region, for simplicity
smooth_sigma = 1 # unit [pixels]
# gaussian smooth broadening parameter, large value will give inaccurate result

# input section, output parameters
coord_print = False
dist_print = False
screen_ini_print = False
screen_fin_print = False

# input section, atomic coordinates, [Angstrom]
coord_inp = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1],
    [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
    [0.5, 0.5, 1], [0.5, 1, 0.5], [1, 0.5, 0.5]
]
# x, y, z

# input section, supercell options, not implemented yet
supercell = False

cell_parameters = [4, 4, 4, 90, 90, 90]
nx = 0
ny = 0
nz = 0

#-------------------------------------------------------------------------------------------
if pxrd:
    runtype_print = 'Powder X-ray Diffraction Simulation'
else:
    runtype_print = '2D X-ray Diffraction Pattern Simulation'

print('\nINITIALIZE| RUN_TYPE: '+runtype_print+'\n'
     +'INITIALIZE| control parameters:\n'
     +'            X-ray wavelength:                           '+str(wvl_inp)+' (nm)\n'
     +'            distance sample placed:                     '+str(dist_inp)+' (m)\n'
     +'            maximum diffraction angle:                  '+str(angle_max_inp)+' (deg)\n'
     +'            resolution level:                           '+str(resol)+'\n'
     +'            coord_scale:                                '+str(coord_scale)+'\n'
     +'            2D diffraction pattern denoise?             '+str(denoise)+'\n'
     +'            Gaussian smooth on 2D pattern?              ['+str(smooth)+']\n'
     +'            use angle resoluted diffraction amplitude?  ['+str(angle_f)+']\n'
     +'            duplicate cell in dimensions?               not implemented\n\n'
     +'INITIALIZE| output parameters:\n'
     +'            print renormalized coordinates?             ['+str(coord_print)+']\n'
     +'            print distance from atom to screen?         ['+str(dist_print)+']\n'
     +'            print intialized screen matrix?             ['+str(screen_ini_print)+']\n'
     +'            print final screen matrix?                  ['+str(screen_fin_print)+']\n'
)

if print_debug_mode:
    exit()

#------------------------------PARAMETERS CONVERTION---------------------------------------
# convert to SI [m]
wvl = wvl_inp * 1E-9
dist = dist_inp

angle_max = float(angle_max_inp)


#------------------------------VARIABLES INITIALIZATION-------------------------------------
scr_edge = dist*math.tan(angle_max*np.pi/180)
print('INITIALIZE| screen width = '+str(scr_edge)+' m')
scr_size = math.ceil(scr_edge) * math.floor(10**resol) # maximum number of elements per dimension
print('INITIALIZE| A '+str(scr_size)+'x'+str(scr_size)+' screen matrix will be created!')

scr = np.zeros((scr_size, scr_size))
if screen_ini_print:
    print('INITIALIZE| screen matrix initailized:')
    print(scr)

scr_dl = scr_edge/scr_size
print('INITIALIZE| pixels coordinates calculated as: \n'
     +'            ix_scr = (index_x-1)*scr_dl, iy_scr = (index_y-1)*scr_dl\n'
     +'            scr_dl = '+str(scr_dl)+' (m)'
)
# propagation distance calculation

    # coordinates centralization...

cell_a = cell_parameters[0]
cell_b = cell_parameters[1]
cell_c = cell_parameters[2]

cell_alpha = cell_parameters[3] # angle between a and b
cell_beta = cell_parameters[4]  # angle between b and c
cell_gamma = cell_parameters[5] # angle between c and a

natom = len(coord_inp)

x_mean = 0
y_mean = 0
z_mean = 0

for iatom in range(natom):

    x_mean += coord_inp[iatom][0]/natom
    y_mean += coord_inp[iatom][1]/natom
    z_mean += coord_inp[iatom][2]/natom

coord_corr = [x_mean, y_mean, z_mean]

coord_cen = np.zeros((natom, 3))

for iatom in range(natom):
    for idim in range(3):

        coord_cen[iatom][idim] = coord_inp[iatom][idim] - coord_corr[idim]

if coord_print:
    print('INITIALIZE| atomic coordinates centeralization:')
    print(coord_cen)

    # distance calculation, distance matrix: [iatom, ix_scr, iy_scr]
dist2scr = np.zeros((natom, scr_size, scr_size))

if matrix_init_debug:
    exit()

#--------------------------------------DISTANCE--------------------------------------------
print('CALCULATE| distance calculation, that may cost a lot of time...')
for iatom in range(natom):
    for idim in range(3):
        for ix_scr in range(scr_size):
            for iy_scr in range(scr_size):
                dist2scr[iatom][ix_scr][iy_scr] = (coord_cen[iatom][0]*10**(coord_scale-11)-ix_scr*scr_dl)**2
                dist2scr[iatom][ix_scr][iy_scr] += (coord_cen[iatom][1]*10**(coord_scale-11)-iy_scr*scr_dl)**2
                dist2scr[iatom][ix_scr][iy_scr] += (coord_cen[iatom][2]*10**(coord_scale-11)-dist)**2
                dist2scr[iatom][ix_scr][iy_scr] = math.sqrt(dist2scr[iatom][ix_scr][iy_scr])
    if dist_print:
        print('CALCULATE| atom '+str(iatom+1)+' distances to screen elements:')
        print(dist2scr[iatom])

#--------------------------------PREPROCESSING OF DIFFRACTION-------------------------------
print('CALCULATE| X-RAY DIFFRACTION STARTS!')
print('CALCULATE| wavelength = '+str(wvl_inp)+' (nm), convert to k, i.e. wavenumber:\n'
     +'           wavenumber k = '+str(1/wvl)+' (m^-1)'
)

if denoise != 'false':
    print('DENOISE| matrix denoise procedure activated!')
    print('DENOISE| denoise method: '+denoise)
    if pxrd:
        print('DENOISE| therefore Powder X-ray Diffraction is not calculated with 2D pattern, but calculated seperately...')
    if denoise == 'cutmean':
        scr_bkgd = 0
    elif denoise == 'cutoff':
        if denoise_cutoff == 0:
            print(' *** ERROR ***: a positive fraction number is expected!')
            exit()
        
if pxrd:
    powder_scr=np.zeros(math.ceil(scr_size*math.sqrt(2)))

#--------------------------------------DIFFRACTION--------------------------------------------
for ix_scr in range(scr_size):
    for iy_scr in range(scr_size):
        for iatom in range(natom):
            if angle_f:
                pixel = math.cos(2*math.pi*dist2scr[iatom][ix_scr][iy_scr]/wvl)
                deg_pixel = math.acos(abs(coord_cen[iatom][2]-dist)/dist2scr[iatom][ix_scr][iy_scr])*math.pi/180 # [deg]
                amp = 1 + angle_res_func[0]*deg_pixel + angle_res_func[1]*deg_pixel**2 +angle_res_func[2]*deg_pixel**3 + angle_res_func[3]*deg_pixel**4
                scr[ix_scr][iy_scr] += amp*pixel
            else:
                scr[ix_scr][iy_scr] += math.cos(2*math.pi*dist2scr[iatom][ix_scr][iy_scr]/wvl)

            if denoise != 'false':
                if denoise == 'cutmean':
                    scr_bkgd += scr[ix_scr][iy_scr]/scr_size/scr_size/natom

            elif pxrd:
                powder_scr[math.floor(math.sqrt(ix_scr**2+iy_scr**2))] += scr[ix_scr][iy_scr]

#---------------------------------------DENOISE------------------------------------------------
if denoise != 'false':
    for ix_scr in range(scr_size):
        for iy_scr in range(scr_size):
            if denoise == 'cutmean':
                if scr[ix_scr][iy_scr] < scr_bkgd*(1+denoise_cutoff):
                    scr[ix_scr][iy_scr] = scr_bkgd
            if denoise == 'cutoff':
                if scr[ix_scr][iy_scr] < denoise_cutoff*np.amax(scr):
                    scr[ix_scr][iy_scr] = 0
            if pxrd:
                powder_scr[math.floor(math.sqrt(ix_scr**2+iy_scr**2))] += scr[ix_scr][iy_scr]
            
#-------------------------------------GAUSSIAN SMOOTH-----------------------------------------
if smooth:
    scr_smooth = np.zeros_like(scr)
    print('SMOOTH| gaussian smooth activated! It will not cost much time but also depends on cutoff you give.')
    for ix_scr in range(scr_size):
        for iy_scr in range(scr_size):
            # every pixel
            for ix_ix_scr in range(max(0, ix_scr-smooth_cutoff), min(ix_scr+smooth_cutoff, scr_size)):
                for iy_iy_scr in range(max(0, iy_scr-smooth_cutoff), min(iy_scr+smooth_cutoff, scr_size)):
                    scr_smooth[ix_scr][iy_scr] += scr[ix_ix_scr][iy_iy_scr]*math.e**(-((ix_ix_scr-ix_scr)**2+(iy_iy_scr-iy_scr)**2)/2/smooth_sigma**2)
    scr = scr_smooth
    print('SMOOTH| smooth complete!')

#-------------------------------------------PLOT----------------------------------------------
if screen_fin_print:
    print('RESULT| screen:')
    print(scr)

if pxrd:
    # usage of subplot: (row, col, index), create grid row*col, index enumerate as:
    # 1 2...
    # n+1 n+2...
    plt.subplot(122)
    pxrd_deg = np.arange(0, angle_max, step=angle_max/len(powder_scr))
    plt.xlabel('diffraction degree theta/deg')
    plt.ylabel('arbitrary unit')
    plt.plot(pxrd_deg[0:len(powder_scr)], powder_scr)

    plt.subplot(121)
    plt.imshow(scr, cmap=plt.cm.gray)
else:
    plt.imshow(scr, cmap=plt.cm.gray)

print('CALCULATE| done.')
plt.show()
