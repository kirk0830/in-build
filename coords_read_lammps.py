import re

# this program read structure file used in lammps
def lammps_input_coord(Infilename):

    try:
        infile = open(Infilename, 'r', encoding='utf-8')
        print('COORD| lammps coordinate file read-in: '+str(Infilename))
    except FileNotFoundError:
        print('COORD| ***error*** lammps coordinate file error: file is non-exist. quit.')
        exit()

    # wash
    # read till find mass information:
    line = infile.readline()
    words = re.split(' |\n', line)
    
    while words[0] != 'Masses':

        line = infile.readline()
        words = re.split(' |\n', line)
    
    print('COORD| #interact# for lammps coordinate files, you need to specify atom types manually. please input their\n'
         +'       element symbols, sperate them with space, so that can be directly used to write standart xyz file.')

    # print till meet atom coordinates:
    while words[0] != 'Atoms':

        line = infile.readline()
        words = re.split(' |\n', line)
        if len(words) >= 3:
            print('       -> element: '+line[:-1])

    elements = input('COORD| waiting for element symbols input: ').split(' ')
    elements = [element for element in elements if element != '']
    # proceed coordinates:
    coord = []
    iloop = 0
    while line and iloop < 50:

        line = infile.readline()
        words = re.split(' |\n', line)
        if len(words) >= 3:

            words = [word for word in words if word != '']
            line2coord = [elements[int(words[1])-1], float(words[-3]), float(words[-2]), float(words[-1])]
            coord.append(line2coord)

    infile.close()

    return coord

def lammps_output_datFile(Infilename):
    pass
