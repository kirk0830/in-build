import numpy as np

class smooth:

    def __init__(self, dataInp = [0]):
        print(
            '=-----------------------------------------------------------------------=\n'
           +'|                                                                       |\n'
           +'|                   Up to two-dimensional SMOOTH code    :P             |\n'
           +'|                                                                       |\n'
           +'=-----------------------------------------------------------------------=\n'
        )

        self.dataInp = dataInp
        self.ndimInp = np.ndim(dataInp)
        self.shapeInp = np.shape(dataInp)

        self.dataOut = np.zeros_like(dataInp)

        if self.ndimInp > 2:
            print(
                'SMOOTH| ***WARNING*** up to now, only datasets that do not exceed 2-\n'
               +'                      dimensions are supported. So only the first 2 \n'
               +'                      dimension will be smoothed.'
                )
            self.ndimInp = 2
        
    def gau1d(self, sigma, fast = False, fast_cutoff = 0):

        print('SMOOTH| 1d-Gaussian function type smooth is activated.')
        if (fast_cutoff - np.floor(fast_cutoff)):
            print('SMOOTH! ***WARNING*** non-integar cutoff value detected, an non-integar cutoff\n'
                 +'                     is not meaningful. I will transfer it using floor method.')
            fast_cutoff = np.floor(fast_cutoff)

        if fast and (int(fast_cutoff)>0):
            print('SMOOTH| fast smooth method is activated, cutoff = '+str(fast_cutoff))

            for ix in range(self.shapeInp[0]):
                for ix_ix in range(max(0, ix-fast_cutoff), min(self.shapeInp[0], ix+fast_cutoff)):
                    self.dataOut[ix] += self.dataInp[ix_ix] * np.e ** (-(ix-ix_ix)**2/2/sigma**2)
        
        else:

            for ix in range(self.shapeInp[0]):
                for ix_ix in range(self.shapeInp[0]):
                    self.dataOut[ix] += self.dataInp[ix_ix] * np.e ** (-(ix-ix_ix)**2/2/sigma**2)

        return self.dataOut

    def fd1d(self, temp = 1e-6, fast = False, fast_cutoff = 0):

        print('SMOOTH| 1d-Fermi-Dirac function type smooth is activated. Smooth temperature = '+str(temp)+' (K)')
        print('SMOOTH| ***ERROR*** however this method is somehow problematic now, I will fix this soon.')
        print('SMOOTH| quit!')
        exit()
        kb = 8.314/6.02E23
        if (fast_cutoff - np.floor(fast_cutoff)):
            print('SMOOTH! ***WARNING*** non-integar cutoff value detected, an non-integar cutoff\n'
                 +'                     is not meaningful. I will transfer it using floor function.')
            fast_cutoff = np.floor(fast_cutoff)

        if fast and (int(fast_cutoff)>0):
            print('SMOOTH| fast smooth method is activated, cutoff = '+str(fast_cutoff))

            for ix in range(self.shapeInp[0]):
                for ix_ix in range(max(0, ix-fast_cutoff), min(self.shapeInp[0], ix+fast_cutoff)):
                    self.dataOut[ix] += self.dataInp[ix_ix] / (1-np.e ** (-abs(ix-ix_ix))/kb/temp)
        
        else:

            for ix in range(self.shapeInp[0]):
                for ix_ix in range(self.shapeInp[0]):
                    self.dataOut[ix] += self.dataInp[ix_ix] / (1-np.e ** (-abs(ix-ix_ix))/kb/temp)

        return self.dataOut

    def gau2d(self, sigma, fast = False, fast_cutoff = 0):

        print('SMOOTH| 2d-Gaussian function type smooth is activated.')
        if (fast_cutoff - np.floor(fast_cutoff)):
            print('SMOOTH! ***WARNING*** non-integar cutoff value detected, an non-integar cutoff\n'
                 +'                     is not meaningful. I will transfer it using floor function.')
            fast_cutoff = np.floor(fast_cutoff)

        if fast and (int(fast_cutoff)>0):
            print('SMOOTH| fast smooth method is activated, cutoff = '+str(fast_cutoff))

            for ix in range(self.shapeInp[0]):
                for iy in range(self.shapeInp[1]):
                    for ix_ix in range(max(0, ix-fast_cutoff), min(self.shapeInp[0], ix+fast_cutoff)):
                        for iy_iy in range(max(0, iy-fast_cutoff), min(self.shapeInp[1], iy+fast_cutoff)):
                            self.dataOut[ix][iy] += self.dataInp[ix_ix][iy_iy] * np.e ** (-((ix-ix_ix)**2+(iy-iy_iy)**2)/2/sigma**2)
        
        else:

            for ix in range(self.shapeInp[0]):
                for iy in range(self.shapeInp[1]):
                    for ix_ix in range(self.shapeInp[0]):
                        for iy_iy in range(self.shapeInp[1]):
                            self.dataOut[ix][iy] += self.dataInp[ix_ix][iy_iy] * np.e ** (-((ix-ix_ix)**2+(iy-iy_iy)**2)/2/sigma**2)

        return self.dataOut
