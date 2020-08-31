import numpy as np
# compensating to numpy, operation like MATLAB
def matrix_minus(A, B):
    # only up to 2-dimensional is supported
    if np.ndim(A) != np.ndim(B):
        raise TypeError
    result = np.zeros_like(A)
    if np.ndim(A) == 2:
        nrow = np.shape(A)[0]
        ncol = np.shape(A)[1]
        for irow in range(nrow):
            for icol in range(ncol):
                result[irow][icol] = A[irow][icol] - B[irow][icol]
    elif np.ndim(A) == 1:
        nele = np.shape(A)[0]
        for iele in range(nele):
            result[iele] = A[iele] - B[iele]
    elif np.ndim(A) == 0:
        result = A - B
    else:
        print('Sorry, other dimensional matrix is not supported yet.')
        raise TypeError
    return result

def matrix_plus(A, B):
    if np.ndim(A) != np.ndim(B):
        raise TypeError
    result = np.zeros_like(A)
    if np.ndim(A) == 2:
        nrow = np.shape(A)[0]
        ncol = np.shape(A)[1]
        for irow in range(nrow):
            for icol in range(ncol):
                result[irow][icol] = A[irow][icol] + B[irow][icol]
    elif np.ndim(A) == 1:
        nele = np.shape(A)[0]
        for iele in range(nele):
            result[iele] = A[iele] + B[iele]
    elif np.ndim(A) == 0:
        result = A + B
    else:
        print('Sorry, other dimensional matrix is not supported yet.')
        raise TypeError
    return result

def matrix_dot_multiply(A, B):
    if np.ndim(A) != np.ndim(B):
        raise TypeError
    result = np.zeros_like(A)
    if np.ndim(A) == 2:
        nrow = np.shape(A)[0]
        ncol = np.shape(A)[1]
        for irow in range(nrow):
            for icol in range(ncol):
                result[irow][icol] = A[irow][icol] * B[irow][icol]
    elif np.ndim(A) == 1:
        nele = np.shape(A)[0]
        for iele in range(nele):
            result[iele] = A[iele] * B[iele]
    elif np.ndim(A) == 0:
        result = A * B
    else:
        print('Sorry, other dimensional matrix is not supported yet.')
        raise TypeError
    return result

def matrix_dot_division(A, B):
    if np.ndim(A) != np.ndim(B):
        raise TypeError
    result = np.zeros_like(A)
    if np.ndim(A) == 2:
        nrow = np.shape(A)[0]
        ncol = np.shape(A)[1]
        for irow in range(nrow):
            for icol in range(ncol):
                result[irow][icol] = A[irow][icol] / B[irow][icol]
    elif np.ndim(A) == 1:
        nele = np.shape(A)[0]
        for iele in range(nele):
            result[iele] = A[iele] / B[iele]
    elif np.ndim(A) == 0:
        result = A / B
    else:
        print('Sorry, other dimensional matrix is not supported yet.')
        raise TypeError
    return result

def matrix_dot_power(A, power):

    result = np.ones_like(A)
    if np.ndim(A) == 2:
        nrow = np.shape(A)[0]
        ncol = np.shape(A)[1]
        for irow in range(nrow):
            for icol in range(ncol):
                ipow = 1
                while ipow <= power:
                    result[irow][icol] *= A[irow][icol]
                    ipow += 1
    elif np.ndim(A) == 1:
        nele = np.shape(A)[0]
        for iele in range(nele):
            ipow = 1
            while ipow <= power:
                result[iele] *= A[iele]
                ipow += 1
    elif np.ndim(A) == 0:
        ipow = 1
        while ipow <= power:
            result *= A
            ipow += 1
    else:
        print('Sorry, other dimensional matrix is not supported yet.')
        raise TypeError
    return result