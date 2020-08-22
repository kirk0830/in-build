import numpy as np

def memoryUpdate(varlog, var):
    varlog.append(var)
    return varlog[-3:]

def LS_simple(xfuncy, x_ini, dx_ini, dx_rescale, conservatism = True, maxloop = 500, verbosity = 'low'):

    """
    # simple minima searching algorithm\n
    This function is used to find local minima by simple accurate line searching method\n
    # input requirement\n
    xfuncy: a function that should be defined external and, just type-in function name to use it in this function\n
    x_ini: one point for initial guess of variable\n
    dx_ini: initial guess of step. a small value is recommanded for very first optimization\n
    dx_rescale: if set to number larger than 1, an aggressive search will be activated and each step is twice of the latter one. If set as 1, a normal searching. If set to a number smaller than 1 but still positive, this search will finally reach a fix point.\n
    conservatism: True is default and recommanded, in case of unknown non-convergence\n
    maxloop: only needed when 'conservatism = True', the max step to find an inteval that may contains minima\n
    verbosity: 'low', 'medium' and 'high' are available options\n
    # output description\n
    a list that contains two element, the first is lower boundary, the second is upper boundary
    """

    # the first parameter must be a function that can yield function value from input x_ini
    x = x_ini
    dx = dx_ini
    try:
        y = xfuncy(x)
    except TypeError:
        print('LS| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    # variables storage
    xlog = [x, x, x]
    ylog = [y, y, y]
    xmem = xlog
    ymem = ylog
    # loop
    iloop = 0
    loopFlag = True

    if verbosity == 'low':
        print('LS| verbosity has been set to \'low\', no information will be printed out on screen.')
    while loopFlag and iloop <= maxloop:

        x += dx
        y = xfuncy(x)
        if verbosity == 'high':
            print('LS| searching point: x = '+str(x)+', y = '+str(y))
        if y > ymem[-1]:
            if iloop == 0:
                # if go larger, back and restart
                if verbosity == 'high':
                    print('LS| search direction is wrong, redirected.')
                dx *= -1
                x = xmem[-1]
                y = xfuncy(x)
            else:
                # should output, temporarily omitted here
                loopFlag = False
                xmem = memoryUpdate(xlog, x)
                ymem = memoryUpdate(ylog, y)
                xout = [xmem[0], xmem[2]]
                if verbosity == 'medium' or verbosity == 'high':
                    print('LS| interval has been found! x_min = '+str(min(xout))+', x_max = '+str(max(xout))+'.')
                return [min(xout), max(xout)]
                
        else:
            # the situation where y_new < y_old
            if verbosity == 'high':
                print('LS| step '+str(iloop)+', larger step is used to find minimium.')
            dx *= dx_rescale
            # save...
            xmem = memoryUpdate(xlog, x)
            ymem = memoryUpdate(ylog, y)
        iloop += 1
        if conservatism == False:
            iloop = 1



