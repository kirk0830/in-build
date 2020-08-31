import numpy as np

def centercoords(coordsinp, masspower = False, totmass = 1, masslist = []):

    """
    This function is to center all atoms. If you want to center them based on their mass, set masspower as True,
    and total mass in relative atomic mass is required, also, mass list is required.\n
    coordsinp: 2 dimensional list, format like standard xyz file\n
    masspower: bool\n
    totmass: float\n
    masslist: 1 dimensional float list, order should be in accord with coordsinp
    """

    print('COORD| i hope you have convert coordinates to float type, if not, no one knows what will happen...')

    natom = np.shape(coordsinp)[0]
    if masspower:
        if np.shape(masslist)[0] != natom:
            print('COORD| ***error*** number of atoms in atomlist is not consistent with masslist. quit.')
            exit()
        else:
            print('COORD| mass-powered atom centering is activated, useful when calculate rotation interia.')

    x_center = 0
    y_center = 0
    z_center = 0

    coordsout = coordsinp

    for iatom in range(natom):

        mass_power = 1

        if masspower:
            mass_power = masslist[iatom]

        x_center += coordsinp[iatom][-3]/natom * (mass_power/totmass)
        y_center += coordsinp[iatom][-2]/natom * (mass_power/totmass)
        z_center += coordsinp[iatom][-1]/natom * (mass_power/totmass)
    
    for iatom in range(natom):

        coordsout[iatom][-3] -= x_center
        coordsout[iatom][-2] -= y_center
        coordsout[iatom][-1] -= z_center

    return coordsout

# one-atom operations

def rotate_atom_x(coord, angle):

    # rotate coordinates around x-axis, which means, rotate coordinates in yz plane.
    # derivation:
    # original vector: rcosA, rsinA
    # rotated vector: rcos(A+B), rsin(A+B) = rcosAcosB - rsinAsinB, rsinAcosB + rcosAsinB
    # re-write in matrix form:
    # 1   0     0  | x       x                       x
    # 0 cosB -sinB | rcosA = rcosAcosB - rsinAsinB = rcos(A+B)
    # 0 sinB  cosB | rsinA   rsinAcosB + rcosAsinB   rsin(A+B)

    angle = angle/180*np.pi

    operator = [
        [1, 0, 0],
        [0, np.cos(angle), -np.sin(angle)],
        [0, np.sin(angle), np.cos(angle)]
    ]

    return np.matmul(coord, operator)

def rotate_atom_y(coord, angle):

    # rotate coordinates around y-axis, which means, rotate coordinates in xz plane.

    # re-write in matrix form:
    # cosB   0   -sinB |
    #  0     1      0  |
    # sinB   0    cosB |

    angle = angle/180*np.pi

    operator = [
        [np.cos(angle), 0, -np.sin(angle)],
        [0, 1, 0],
        [np.sin(angle), 0, np.cos(angle)]
    ]

    return np.matmul(coord, operator)

def inversion_atom(coord, mode = 'xyz'):

    if mode == 'xyz':
        coord[-3] = -coord[-3] # x
        coord[-2] = -coord[-2] # y
        coord[-1] = -coord[-1] # z
    elif mode == 'x':
        # same as reflection towards yz-plane.
        coord[-3] = -coord[-3]
    elif mode == 'y':
        # same as reflection towards xz-plane.
        coord[-2] = -coord[-2]
    elif mode == 'z':
        # same as reflection towards xy-plane.
        coord[-1] = -coord[-1]
    elif mode == 'xy':
        coord[-3] = -coord[-3] # x
        coord[-2] = -coord[-2] # y
    elif mode == 'yz':
        coord[-2] = -coord[-2] # y
        coord[-1] = -coord[-1] # z
    elif mode == 'xz':
        coord[-3] = -coord[-3] # x
        coord[-1] = -coord[-1] # z
    else:
        exit()
    return coord

def mirror_xy_atom(coord):

    return inversion_atom(coord, mode='z')

def mirror_xz_atom(coord):

    return inversion_atom(coord, mode='y')

def mirror_yz_atom(coord):

    return inversion_atom(coord, mode='x')

# poly-atom operations
# visit the 1st atom in list1, find atom_i in list2, residue: list2[0:i-1], list2[i+1:], add together name as list2
# visit the 2nd atom in list1, find atom_j in list2...
# ...

def findinlist(atom, atomlist):

    # do not call this subroutine from external programs except you know what you are doing.

    natom = np.shape(atomlist)[0]
    nTrue = 0
    index = -1

    for iatom in range(natom):
        if np.array_equal(atom, atomlist[iatom][:]):
            nTrue += 1
            index = iatom

    if nTrue > 1:
        print('COORD| ***error*** atoms overlap! quit.')
        exit()
    else:
        return index

def isCoordSame(coords1, coords2, acc = 6):

    natom = np.shape(coords1)[0]

    if np.shape(coords2)[0] != natom:
        print('COORD| ***error*** atom number is not in consistency.')
        return False
    else:
        print('COORD| atom number check passed. :) Good luck.')
    # pre-processing of coordinates
    for iatom in range(natom):

        coords1[iatom][-1] = round(coords1[iatom][-1], acc)
        coords1[iatom][-2] = round(coords1[iatom][-2], acc)
        coords1[iatom][-3] = round(coords1[iatom][-3], acc)
        coords2[iatom][-1] = round(coords2[iatom][-1], acc)
        coords2[iatom][-2] = round(coords2[iatom][-2], acc)
        coords2[iatom][-3] = round(coords2[iatom][-3], acc)
    
    print('COORD| convert two atomlists to the same precision level: '+str(10**-acc))
    list2find = coords2
    index = -2

    identicalFlag = False
    for iatom in range(natom):

        if index == -2:
            index = findinlist(coords1[iatom][:], list2find)
                
        if iatom > 0:
            if index >= 0:
                list2find = list2find[0:index]+list2find[index+1:]
                index = findinlist(coords1[iatom][:], list2find)
                if iatom == natom-1 and index >= 0:
                    identicalFlag = True
            else:
                break

    return identicalFlag

def cell_duplication(coords, a, b, c, alpha = 90, beta = 90, gamma = 90, nx=1, ny=1, nz=1):
    print('SUPERCELL| supercell function is activated. I hope you have type-in parameters\n'
         +'           that all required. alpha, beta, gamma will use 90 deg as default if \n'
         +'           not explicitly defined.'
    )
    if np.ndim(coords) != 2:
        print('SUPERCELL| ***ERROR*** may be you forget to give coordinates list. Program is terminated.')
        exit()

    if a<=0 or b<=0 or c<=0 or alpha<=0 or beta<=0 or gamma<=0 or alpha>=180 or beta>=180 or gamma>=180:
        print('SUPERCELL| ***ERROR*** invalid parameter input! Program is terminated.')
        exit()
    
    a = float(a)
    b = float(b)
    c = float(c)
    alpha_r = float(alpha/180) * np.pi
    beta_r = float(beta/180) * np.pi
    gamma_r = float(gamma/180) * np.pi
    a_vec = [a, 0, 0]
    b_vec = [b*np.cos(alpha_r), b*np.sin(alpha_r), 0]
    c1 = c*np.cos(gamma_r)
    c2 = c*(np.cos(beta_r)-np.cos(alpha_r)*np.cos(gamma_r))/np.sin(alpha_r)
    c3 = np.sqrt(c**2 - c1**2 - c2**2)
    c_vec = [c1, c2, c3]

    print('SUPERCELL| cell parameters information:\n'
         +'           a = '+str(a)+' Angstrom | alpha = '+str(alpha)+' deg\n'
         +'           b = '+str(b)+' Angstrom | beta = '+str(beta)+' deg\n'
         +'           c = '+str(c)+' Angstrom | gamma = '+str(gamma)+' deg\n'
         +'           cell(1) = '+str(a_vec)+'\n'
         +'           cell(2) = '+str(b_vec)+'\n'
         +'           cell(3) = '+str(c_vec)+'\n'
    )
    
    natom = np.shape(coords)[0]
    print('SUPERCELL| original cell information: '+str(natom)+' atoms in total.')
    ncol = np.shape(coords)[1]

    if nx != int(nx) or ny != int(ny) or nz != int(nz) or nx <= 0 or ny <= 0 or nz <= 0:
        print('SUPERCELL| invalid duplicate number, they should have been integars. quit.')
        exit()
    
    coords_readin = coords
    coords_add = []
    if ncol == 3:
        print('SUPERCELL| there seems no element information in coordinates input, single-\n'
             +'           element treatment will be employed...'
        )
        for iatom in range(natom):
            for inx in range(nx):
                for iny in range(ny):
                    for inz in range(nz):
                        atom_dup = coords[iatom][:]
                        atom_dup[0] += (inx*a_vec[0] + iny*b_vec[0] + inz*c_vec[0])
                        atom_dup[1] += (inx*a_vec[1] + iny*b_vec[1] + inz*c_vec[1])
                        atom_dup[2] += (inx*a_vec[2] + iny*b_vec[2] + inz*c_vec[2])
                        if inx or iny or inz:
                            coords_add.append(atom_dup)
    elif ncol == 4:
        for iatom in range(natom):
            for inx in range(nx):
                for iny in range(ny):
                    for inz in range(nz):
                        # atom_dup = coords[iatom] # do not execute this line!
                        atom_dup = coords[iatom][:]
                        atom_dup[1] += (inx*a_vec[0] + iny*b_vec[0] + inz*c_vec[0])
                        atom_dup[2] += (inx*a_vec[1] + iny*b_vec[1] + inz*c_vec[1])
                        atom_dup[3] += (inx*a_vec[2] + iny*b_vec[2] + inz*c_vec[2])
                        if inx or iny or inz:
                            coords_add.append(atom_dup)
# developer notes here: it is interesting that if use:
# atom_dup = coords[iatom], atom_dup is not a new object, instead, it is, the constant pointer,
# or say, the re-named coords[iatom]. So latter steps will directly change data in coords[][],
# which is not what we expect.
# instead, if we write atom_dup = coords[iatom][:], a new object will be created, all following
# steps will proceed as usual.
    else:
        print('SUPERCELL| ***ERROR*** wrong coordinates format. Program is terminated.')
        exit()

    return coords_readin+coords_add

def coords_rotate(
    coords, 
    angleX = 0, 
    angleY = 0
    ):

    natom = np.shape(coords)[0]
    for iatom in range(natom):
        coord = coords[iatom][1:]
        coord = rotate_atom_x(coord, angleX)
        coord = rotate_atom_y(coord, angleY)
        coords[iatom][-3] = coord[0]
        coords[iatom][-2] = coord[1]
        coords[iatom][-1] = coord[2]
    return coords

def coords_sort_prep(coordlist):

    rows = np.shape(coordlist)[0]

    coordlist_prep = []
    for irow in range(rows):

        coordlist_prep.append(
            (
                coordlist[irow][0], 
                coordlist[irow][1], 
                coordlist[irow][2], 
                coordlist[irow][3]
                )
                )
    return coordlist_prep

def coords_sort_postp(coordlist_sorted):

    rows = np.shape(coordlist_sorted)[0]

    coordlist_postp = []
    for irow in range(rows):

        coordlist_postp.append(
            [
                str(coordlist_sorted[irow][0]).split('\'')[1], 
                coordlist_sorted[irow][1], 
                coordlist_sorted[irow][2], 
                coordlist_sorted[irow][3]
                ]
                )

    return coordlist_postp

def coords_sort(
    coordlist, 
    axis, 
    multisort = False, 
    axislist = []
    ):

    header = [('element', 'S10'), ('x', float), ('y', float), ('z', float)]

    print('COORD| coordinates sorting operation is activated.\n'
         +'       ***warning*** three-column coordinates data is not supported, please, add dumb elements before using.')
    # developer note:
    # there is a bug that, numpy_bytes is obtained if directly use commands stated in the following. Which will add a
    # b'str' at the head of original string. So one must convert explicitly it into str data type and extract proper
    # contents from it, like what will do in function coords_sort_postp().
    coordlist_prep = coords_sort_prep(coordlist)
    entitled = np.array(coordlist_prep, dtype = header)


    if axis == 'x' or axis == '1':

        coord_sorted = np.sort(entitled, order='x')
    elif axis == 'y' or axis == '2':

        coord_sorted = np.sort(entitled, order='y')
    elif axis == 'z' or axis == '3':

        coord_sorted = np.sort(entitled, order='z')
    else:

        print('COORD| ***error*** axis parameter can not be recognized, please check, quit.')
        exit()

    coordlist_postp = coords_sort_postp(coord_sorted)

    return coordlist_postp

# batch substitution
def coords_batchSubs_plane(
    fromElement, 
    toElement, 
    coordlist, 
    axis, 
    low, 
    high
    ):

    """
    # Input requirement\n
    fromElement: element that choose to be substituted, if set to 'all', all atoms in region will be substituted\n
    toElement: element that choose to overwrite\n
    coordlist: original coordinate list, 2-dimensional list, 4 columns, str, float, float, float type respectively\n
    axis: 'x', 'y' or 'z', in future, will be modified as 'a', 'b' or 'c' that compatible with all shaped cell\n
    low: lower boundary of region want to substitute\n
    high: upper boundary of region want to substitute
    """
    
    rows = np.shape(coordlist)[0]
    if axis == 'x':
        axis = 1
    elif axis == 'y':
        axis = 2
    elif axis == 'z':
        axis = 3
    elif axis == 'a':
        print('COORD| ***error*** presently non-orthdo. cell is not supported yet, quit.')
        exit()
    elif axis == 'b':
        print('COORD| ***error*** presently non-orthdo. cell is not supported yet, quit.')
        exit()
    elif axis == 'c':
        print('COORD| ***error*** presently non-orthdo. cell is not supported yet, quit.')
        exit()
    else:
        print('COORD| ***error*** bad axis requirement, please check your input, quit.')
        exit()
    
    for irow in range(rows):

        if (
            (fromElement == 'all' or fromElement == coordlist[irow][0])
             and 
             (coordlist[irow][axis] >= low and coordlist[irow][axis] <= high)
             ):
            coordlist[irow][0] = toElement

    return coordlist

def coords_batchSubs_sphere(
    fromElement, 
    toElement, 
    coordlist, 
    center = [0, 0, 0], 
    radius = -1,
    shell = False,
    thickness = -1
    ):
    """
    # Input requirement\n
    fromElement: element that choose to be substituted, if set to 'all', all atoms in region will be substituted\n
    toElement: element that choose to overwrite\n
    coordlist: original coordinate list, 2-dimensional list, 4 columns, str, float, float, float type respectively\n
    center: sphere center coordinates, 1-dimensional list is expected\n
    radius: radius of spherical region, if shell mode is enabled, this is outer radius\n
    shell: if set to true, thickness must be specified, atoms residing between radius and radius+thickness will be substituted\n
    thickness: thickness of spherical shell, must be positive
    """

    if radius <= -1:
        print('COORD| ***error*** unreasonable radius value detected, radius must be specified with a positive value. quit.')
        exit()
    if shell == False:
        thickness = radius
    radius1 = radius - thickness
    radius2 = radius

    rows = np.shape(coordlist)[0]

    for irow in range(rows):

        dist = np.sqrt(
            (center[0]-coordlist[irow][1])**2
           +(center[1]-coordlist[irow][2])**2
           +(center[2]-coordlist[irow][3])**2
        )
        if (
            (fromElement == 'all' or fromElement == coordlist[irow][0])
            and
            (dist >= radius1 and dist <= radius2)
        ):
            coordlist[irow][0] = toElement
    
    return coordlist

def coords_batchSubs_cube(
    fromElement, 
    toElement, 
    coordlist, 
    point1, 
    point2,
    shell = False,
    point3 = [],
    point4 = []
    ):
    

    pass

def coords_truncate_sphere(
    coordlist, 
    center = [0, 0, 0], 
    radius = 0, 
    domain = 'out', 
    verbosity = 'low'
    ):
    # to remove all atoms in selected region
    natom = np.shape(coordlist)[0]
    if verbosity == 'high' or verbosity == 'debug':
        print('COORD| coords_transform_truncation operation activated, atom counting: '+str(natom))
    # pbc not supported yet
    # developer note: again, here, one should not directly write "coordsout = coordlist", or at latter steps,
    # operation on coordsout will directly copy to coordlist, for again, coordsout is just a rename of coordlist.
    coordsout = coordlist[:][:]
    for iatom in range(natom):
        try:
            if verbosity == 'debug':
                print('COORD| distance calculation for truncation, atom '+str(iatom)+' coordinate: '+str(coordlist[iatom][1:]))
            dist = np.sqrt(
                (coordlist[iatom][-3] - center[0])**2
               +(coordlist[iatom][-2] - center[1])**2
               +(coordlist[iatom][-1] - center[2])**2
                )
            if (domain == 'out' and dist > radius) or (domain == 'in' and dist < radius):
                if verbosity == 'high' or verbosity == 'debug':
                    print('COORD| remove atom '+str(iatom))
                coordsout.remove(coordlist[iatom][:])
        except IndexError:
            print('COORD| ***warning*** please check if anything wrong with data of atom '+str(iatom))
    
    if verbosity == 'high' or verbosity == 'debug':
        print('COORD| truncation complete, number of atoms left: '+str(np.shape(coordsout)[0]))
    return coordsout

