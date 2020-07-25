import numpy as np

def dup(coords, a, b, c, alpha = 90, beta = 90, gamma = 90, nx=1, ny=1, nz=1):
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