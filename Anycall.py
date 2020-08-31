from coords_io import readcoords, writecoords
from coords_transform import coords_batchSubs_plane
import numpy as np

xyz = readcoords(filename = 'Pd864-202008311546.xyz', convert = True)

xyzO = coords_batchSubs_plane(fromElement = 'Pd', toElement = 'O', coordlist = xyz, axis = 'x', low = 0, high = 4)
writecoords(coordlist = xyzO, projectname = 'Pd_plane')