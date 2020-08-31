from distance import dist
import numpy as np

# this function is to generate distance diagnoal matrix
def distTriangleMat(
                    coordslist, 
                    full_matrix = True,
                    pbc = False,
                    cell_a = 0,
                    cell_b = 0,
                    cell_c = 0,
                    cell_alpha = 90,
                    cell_beta = 90,
                    cell_gamma = 90
                    ):
    
    """
    function to calculate distance between ALL atoms, return a OD matrix, 2 dimensional.\n
    # input requirement\n
    coordlist: atomic coordinates, standardized and widely used in this package\n
    full_matrix: [bool] output a fully diagonal matrix or not. If set as false, half of elements in output matrix will be zero cuz have not been calculated\n
    pbc: [bool] if use periodic boundary condition correction. Useful in periodic system, set to False if called in an isolated system calculation\n
    other six parameters: MUST be set reasonable values if set pbc as True, these are cell parameters\n
    """

    natom = np.shape(coordslist)[0]
    # dont know how many columns the coordlist has, either 3 or 4, use -1, -2 and -3 to visit
    distmat = np.zeros((natom, natom))
    print('DISTANCE| distance matrix generated, size: '+str(natom)+'x'+str(natom)+'.')
    for i in range(natom):

        if pbc:
            distlist = dist(
                            atomlist=coordslist[i+1:][:], 
                            pointxyz=coordslist[i][-3:],
                            point_in_coord=True,
                            pbcFlag=True,
                            pbc_a=cell_a,
                            pbc_b=cell_b,
                            pbc_c=cell_c,
                            pbc_alpha=cell_alpha,
                            pbc_beta=cell_beta,
                            pbc_gamma=cell_gamma
                            )
        
        else:
            distlist = dist(
                            atomlist=coordslist[i+1:][:], 
                            pointxyz=coordslist[i][-3:]
                            )

        print('DISTANCE| calculate distance list from atom '+str(i)+' to other atoms...')

        for j in range(natom-i-1):

            distmat[i][i+j+1] = distlist[j]
            if full_matrix:
                distmat[i+j+1][i] = distmat[i][i+j+1]

    return distmat

# for a given cutoff value, return the dict that contains nearest neighboring particles' indices
def nearestNeighbors(distmat, cutoff):

    print('PAIR| count nearest neighboring particles... current cutoff value in Angstrom: '+str(cutoff))
    natom = np.shape(distmat)[0]

    # initialization of neighborDict:
    neighbordict = {}
    for i in range(natom):

        neighbordict[i] = []

    for i in range(natom):
        for j in range(natom-i-1):

            if distmat[i][i+j+1] < cutoff:

                # structured data
                neighbordict[i].append((i+j+1, distmat[i][i+j+1]))
                neighbordict[i+j+1].append((i, distmat[i][i+j+1]))

    print('PAIR| nearestNeighbors information has been saved in dict. Note that data in parenthesis has\n'
         +'      no difference with that in bracket. It is just for convienence of later pair-calculation.')
    return neighbordict

def triangleExist(pointA, pointB, pointC, collection):

    aTriangle = [pointA, pointB, pointC]
    nComb = np.shape(collection)[0]

    for iComb in range(nComb):

        if np.array_equal(aTriangle, collection[iComb]):

            return True

def q_tetra_order_parameter(coords, cutoff):

    """
    q_tetra order parameter calculation\n
    # Definition\n
    q_tetra of atom i is denoted as q_tetra_i,\n
    q_tetra_i = 1 - 3/8 sum((cos(theta)+1/3)^2, theta), theta)\n
    theta is the angle of j-i-k, j and k are atoms neighboring to atom i, the summation above goes over all possible atoms j and k.
    # Input requirement\n
    coords: atomic coordinates, 2 dimensional list, as used over all functions in this package\n
    cutoff: cutoff value given manually to calculate only neigboring atoms within a certain distance, unit is Angstrom\n
    # Output description\n
    a dictionary that contains all possible atoms' q_tetra
    """
    print('PAIR| q_tetra order parameter analysis is activated, a series of calculations will be carried out.')
    distMatrix = distTriangleMat(coordslist=coords)
    neighors = nearestNeighbors(distmat=distMatrix, cutoff=cutoff)

    description = [('index',int),('distance',float)]
    natom = len(neighors)

    q_tetra_dict = {}
    for iatom in range(natom):

        connectAtoms = neighors[iatom]
        q_tetra_thisAtom = 0

        # format: [(1, dist1), (2, dist2), ...]
        if len(connectAtoms) >= 4:
            print('PAIR| q_tetra-calculation is performed on atom '+str(iatom))
            data_to_sort = np.array(connectAtoms, dtype=description)
            data_sorted = np.sort(data_to_sort, order='distance')
            points_to_calc = [data_sorted[0][0], data_sorted[1][0], data_sorted[2][0], data_sorted[3][0]]
            print('PAIR| q_tetra: atom '+str(iatom)+' -> atoms '+str(points_to_calc))
            # start
            for atom_i in range(3):
                for disp_i in range(1,4-atom_i):
                    # vector_i
                    vector_i = [
                        coords[points_to_calc[atom_i]][-3]-coords[iatom][-3],
                        coords[points_to_calc[atom_i]][-2]-coords[iatom][-2],
                        coords[points_to_calc[atom_i]][-1]-coords[iatom][-1]
                        ]
                    vector_j = [
                        coords[points_to_calc[atom_i+disp_i]][-3]-coords[iatom][-3],
                        coords[points_to_calc[atom_i+disp_i]][-2]-coords[iatom][-2],
                        coords[points_to_calc[atom_i+disp_i]][-1]-coords[iatom][-1]  
                    ]
                    cos_ij = np.dot(vector_i, vector_j)/np.linalg.norm(vector_i)/np.linalg.norm(vector_j)
                    q_tetra_thisAtom -= 3/8*(cos_ij+1/3)**2

            q_tetra_thisAtom += 1
            q_tetra_dict[iatom] = q_tetra_thisAtom

        else:
            print('PAIR| for atom '+str(iatom)+' has neighboring atoms less than 4, q_tetra-calculation skips.')

    return q_tetra_dict

# ------------------------------more complicated version of order parameters------------------------------
from orthopoly.spherical_harmonic import sph_har
#from orthopoly.spherical_harmonic import cart2sph
# orthopoly.spherical_harmonic.cart2sph() has fatal bug: I can not understand why vector (1, 0, 0)
# will be converted into (1, pi/2, pi), which means theta = 0 is from -x, it is, not a common definition
# so I write by myself instead ->

def cart2sph(x, y, z):

    r = np.sqrt(x**2+y**2+z**2)
    if r <= 0:
        print('PAIR| ***error*** bad distance calculation!!! distance is an non-positive value, quit.')
        exit()
    # zeroDivision and negative value exclusion
    # [r, np.arccos(z/r), np.arctan(y/x)]
    if x*y < 0:
        return [r, np.arccos(z/r), 2*np.pi+np.arctan(y/x)]
    if x*y > 0:
        return [r, np.arccos(z/r), np.arctan(y/x)]
    if x*y == 0:
        if x == 0:
            if y > 0:
                return [r, np.arccos(z/r), np.pi/2]
            if y < 0:
                return [r, np.arccos(z/r), 3*np.pi/2]
            if y == 0:
                return [r, 0, 0]
        if x > 0:
            return [r, np.arccos(z/r), 0]
        if x < 0:
            return [r, np.arccos(z/r), np.pi]
        
# this package can be safely used.
# formulation reference:
# English version and package document: 
# https://wordsworthgroup.github.io/orthopoly/spherical_harmonic.html
# Chinese version of spherical harmonics chart: 
# https://baike.baidu.com/item/%E7%90%83%E8%B0%90%E5%87%BD%E6%95%B0%E8%A1%A8/22781995?fr=aladdin

# developer note: although this is somewhat slow and I do not know why.

def q_l_order_parameter(
    coords, 
    cutoff, 
    l = 6, 
    rot_invariant = False, 
    pbc = False, 
    verbosity = 'low'):

    '''
    # q_l order parameter calculation\n
    To calculate local structure that has special rotation symmetry, as first demonstrated in 10.1103/PhysRevB.28.784\n
    Formulation refers to https://www.nature.com/articles/srep00505.pdf\n
    l = 3, tetrahedron order -> q_tetra is independently implemented previously.\n
    l = 6, close-packed, fcc or hcp\n
    l = 8, cubic order\n
    # formulation\n
    coordinates-dependent: q_lm(i) = 1/N_neibghor * sum{over all neighbors, denote as j}(Y_lm(r_ij))\n
    note: q_lm(i) is an element of (2l+1)-dimensional vector, where m is the index from -l to l.\n
    rotation-invariant: q_l(i) = sqrt(4*pi/(2*l+1)*norm(q_lm(i)))\n
    # Input requirement\n
    coords: coordinates list in 2-dimensional that universaly used in SIMUPKGS\n
    l: angular momentum number, determining symmetry of spherical harmonics\n
    rot_invariant: calculate norm of q\n
    pbc: periodic boundary condition correction, if enabled, cutoff value will be used to generate images of atoms
    residing far away from presently examined atom. This tag will be implemented in the future\n
    verbosity: 'low', 'medium', 'high' and 'debug' print level is available.
    '''
    print('PAIR| q_l order parmeter analysis is activated, present selection: q'+str(l))

    # symmetry property parameters
    if l == 6:
        # fcc or hcp
        N_neighbor = 12
    elif l == 3:
        # tetrahedron
        N_neighbor = 4

    # nearest neighbors selection...
    distMatrix = distTriangleMat(coordslist=coords)
    neighors = nearestNeighbors(distmat=distMatrix, cutoff=cutoff)
    # header of structral list
    description = [('index',int),('distance',float)]
    natom = len(neighors)

    ndim = 2*l+1
    q = np.zeros((natom, ndim))

    for iatom in range(natom):
        # for all atoms in coordinate file...
        # neighboring list without sorting, in mess
        connectAtoms = neighors[iatom]

        # format: [(1, dist1), (2, dist2), ...]
        if len(connectAtoms) >= N_neighbor:

            # the atom that indeed has over than N_neighor neighboring atoms
            # if N_neighbor is a large number, one should use large cutoff value to have a completed neighborlist
            if verbosity == 'high':
                print('PAIR| q-order parameter calculation is performed on atom '+str(iatom))

            data_to_sort = np.array(connectAtoms, dtype=description)
            if verbosity == 'debug':
                print('Atom '+str(iatom)+' data_to_sort: '+str(data_to_sort))
            data_sorted = np.sort(data_to_sort, order='distance')
            if verbosity == 'debug':
                print('Atom '+str(iatom)+' data_sorted: '+str(data_sorted))
            points_to_calc = []
            for ineighbor in range(N_neighbor):
                points_to_calc.append(data_sorted[ineighbor][0])
            if verbosity != 'low':
                print('PAIR| q_l: atom '+str(iatom)+' -> atoms '+str(points_to_calc))
            # nearest-neighboring atoms saved, holding atomic indices in points_to_calc list

            # start
            if verbosity != 'low':
                print('PAIR| (2l+1)-dimensional q-vector is calculated...')
            
            for idim in range(ndim):
                for iNeighbor in range(N_neighbor):
                    r = [
                        coords[points_to_calc[iNeighbor]][-3] - coords[iatom][-3],
                        coords[points_to_calc[iNeighbor]][-2] - coords[iatom][-2],
                        coords[points_to_calc[iNeighbor]][-1] - coords[iatom][-1]
                    ]
                    polarAngle = cart2sph(r[0], r[1], r[2])[1:]
                    if verbosity == 'debug':
                        print('PAIR| cart2polar conversion, polar angles in rad: (atom '+str(iatom)+' )'+str(polarAngle))
                    q[iatom][idim] += 1/N_neighbor*sph_har(t = polarAngle[0], p = polarAngle[1], n = l, m = idim-l)
                    if verbosity == 'debug':
                        print('PAIR| single q-value calculation passed...')
        else:
            if verbosity != 'low':
                print('PAIR| for atom '+str(iatom)+' has neighboring atoms less than required, q-calculation skips.')

    if rot_invariant:
        if verbosity != 'low':
            print('PAIR| rotation-invariant number is required...')
        q_invar = np.zeros(natom)
        for iatom in range(natom):
            q_invar[iatom] = np.sqrt((4*np.pi/(2*l+1))*np.linalg.norm(q[iatom][:]))
        return [q, q_invar]
    else:
        return q
