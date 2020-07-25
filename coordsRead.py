import numpy as np

def read(filename = None, convert = False):

    coords_inp = []
    print('COORD| coordinates read-in function is activated. Parsing file: '+filename)
    if filename == None:
        print('COORD| ***ERROR*** no coordinate file is read, program is terminated.')
        exit()

    fileTag = open(filename, 'r', encoding='utf-8')
    print('COORD| the first two lines in *.xyz file will be omitted. Although it may be a\n'
         +'       good way to check completeness of file.')

    lineskip = 2
    istep = 1
    while istep <= lineskip:
        iline = fileTag.readline()
        istep += 1

    while(len(iline)>3):
        iline = fileTag.readline()[:-1].split(' ')
        iline = [iwords for iwords in iline if iwords != '']
        if len(iline):
            coords_inp.append(iline)

    fileTag.close()

    if convert:
        print('COORD| read-in complete, converting coordinates from str to float. Here atoms will be counted.')
        coordShape = np.shape(coords_inp)
        print('COORD| there are '+str(coordShape[0])+' atoms in total.')

        for iatom in range(coordShape[0]):
            for iaxis in range(1,4):
                coords_inp[iatom][iaxis] = float(coords_inp[iatom][iaxis])
        
        print('COORD| conversion complete.')
    
    else:
        print('COORD| read-in complete, remember that convert data from str to float-type for future use.')

    return coords_inp
