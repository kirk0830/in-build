# Developer note: ALL functions created below will be extended to support multi-dimensional function
# optimization, up to now, only 1-dimensional functions are supported.

import numpy as np
import _in_matrix_op as mop

def memoryUpdate(varlog, var):
    varlog.append(var)
    return varlog[-3:]

def direct_search(
    xfuncy, 
    x_ini, 
    dx_ini, 
    dx_rescale, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):

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

        x = mop.matrix_plus(x, dx)
        y = xfuncy(x)
        if verbosity == 'high':
            print('LS| searching point: x = '+str(x)+', y = '+str(y))
        if y > ymem[-1]:
            if iloop == 0:
                # if go larger, back and restart
                if verbosity == 'high':
                    print('LS| search direction is wrong, redirected.')
                dx = np.dot(-1, dx)
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
            dx = np.dot(dx_rescale, dx)
            # save...
            xmem = memoryUpdate(xlog, x)
            ymem = memoryUpdate(ylog, y)
        iloop += 1
        if conservatism == False:
            iloop = 1

def golden_search(
    xfuncy, 
    x_min, 
    x_max, 
    epsilon, 
    conservatism = True, 
    maxloop = 500, 
    otherRatio = False, 
    ratio = 0.618, 
    verbosity = 'low'
    ):

    """
    # Golden search algorithm\n
    a simple method that finds one local minima.\n
    # Input requirement\n
    xfuncy: external function of y(x)\n
    x_min: lower boundary of inteval where to find minima\n
    x_max: upper boundary of inteval where to find minima\n
    epsilon: accurancy control parameter, the larger you set, the more precise you will get\n
    conservatism: see document of function "direct_search" function\n
    maxloop: see document of function "direct_search" function\n
    otherRatio: [bool] although this function is named Golden_search, take in mind that Golden just means 0.618, so other proportion is also supported\n
    ratio: [float] if set otherRatio = True, you can use an arbitrary value to shrink your inteval\n
    # Output description\n
    a list contains both modified lower and upper boundaries
    """
    x_left = x_min
    x_right = x_max
    try:
        y_left = xfuncy(x_left)
        y_right = xfuncy(x_right)
    except TypeError:
        print('LS| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    
    if verbosity != 'low':
        print('LS| golden search method is activated, search domain: ['+str(x_min)+', '+str(x_max)+'].')
    
    x_left_mem = []
    x_right_mem = []
    y_left_mem = []
    y_right_mem = []

    if otherRatio:
        t = ratio
        if verbosity != 'low':
            print('LS| golden echo: a modified ratio to segment inteval is used. ratio input: '+str(t))
    else:
        t = 0.618
    
    # we use epsilon as the magnitude of inteval, i.e., if epsilon = 1, convergence requirement: x_max - x_min <= 1E-1*(x_max_0 - x_min_0)
    iloop = 0
    while (iloop <= maxloop) and (
        np.linalg.norm(mop.matrix_minus(x_right, x_left))
        >= 
        10**(-epsilon)*np.linalg.norm(mop.matrix_minus(x_max, x_min))
        ):
        iloop += 1

        if y_right >= y_left:
            x_right_mem.append(x_right)
            x_right = mop.matrix_plus(
                x_left, 
                np.dot(max(t, 1-t), mop.matrix_minus(x_right, x_left))
                )
            y_right = xfuncy(x_right)
            y_right_mem.append(y_right)
            if verbosity == 'high':
                print('LS| shrink inteval upper boundary: '+str(x_right_mem[-1])+' -> '+str(x_right))
        else:
            x_left_mem.append(x_left)
            x_left = mop.matrix_plus(
                x_left, 
                np.dot(min(t, 1-t), mop.matrix_minus(x_right, x_left))
                )
            y_left = xfuncy(x_left)
            y_left_mem.append(y_left)
            if verbosity == 'high':
                print('LS| shrink inteval lower boundary: '+str(x_left_mem[-1])+' -> '+str(x_left))
        
        if verbosity != 'low':
            print('LS| present inteval: ['+str(x_left)+', '+str(x_right)+'].')
        if conservatism == False:
            iloop = 1

    return [x_left, x_right]

def quadratic_search(
    xfuncy, 
    x_ini, 
    dx_ini, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):
    """
    # Usage Warning\n
    This function is based on consequential 2-spline fitting. So it may always happen if local curvature is negative that this function will find a maximum instead of
    minimum, so this function is not recommended and will only coded in the future. Instead, I will most probably write a ploynomial_search that supports 2- and 3-spline
    both.
    """
    pass

def armijo_fuzzy_search_1d(
    xfuncy, 
    xGrady, 
    x_ini, 
    dx_ini, 
    dxInterface = False,
    sigma = 0.2, 
    beta = 0.5, 
    converLevel = 20, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):
    """
    # Original Armijo fuzzy search algorithm\n
    This is original version of Armijo algorithm, where step in it is quite arbitrary and may be unreasonable under some circumstances.\n
    # Formulation (Armijo criterion)\n
    Given parameters beta in (0, 1) and sigma in (0, 0.5), for m belonging to non-negative integars:\n
    if inequality f(x_k + beta**m * dx) <= f(x_k) + sigma * beta**m * Df(x_k) * dx can be satisfied,\n
    use the minimal m as output and update x_k <- x_k + beta**m * dx is admitted.\n
    This criterion is equivalent with:\n
    [f(x_k + beta**m * dx) - f(x_k)]/(beta**m * Df(x_k) * dx) ~=~ Df(x_k) <= sigma*Df(x_k),\n
    that is, ensuring Df(x_k) is negative.\n
    # Input requirement\n
    xfuncy: external function y(x)\n
    xGrady: derivative of y(x)\n
    x_ini: starting point of x\n
    dx_ini: for internal usage, i.e., dxInterface = False, updating step magnitude, dx = 1, stepsize = 1E-1*x_ini;
    for external usage that dx provided by other algorithm, set dxInterface = True and dx is directly given the proper value.\n
    sigma: derivative shrinking coefficent\n
    beta: step shrinking coeffient\n
    converLevel: maximum of order of m to try.
    """

    # Type chcek:
    try:
        y_ini = xfuncy(x_ini)
        yy_ini = xGrady(x_ini)
        if verbosity == 'high':
            print('LS| initial point information: x = '+str(x_ini)+', y = '+str(y_ini)+', dy/dx = '+str(yy_ini))
    except TypeError:
        print('LS| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    
    if dxInterface:
        dx = dx_ini
    else:
        dx = x_ini * 10**(-dx_ini) # will be deprecated in advanced oprimization methods
    
    x = x_ini
    x_mem = [x]
    m_mem = [] 

    iloop = 0
    if_reversed_dx = False
    while iloop <= maxloop:

        m=0

        while xfuncy(x+beta**m*dx) >= (xfuncy(x) + sigma*(beta**m*dx)*xGrady(x)) and m<=converLevel:
            print('LS| Armijo order m = '+str(m)
                 +', terms of inequality: '+str(xfuncy(x+beta**m*dx))+', '+str(xfuncy(x) + sigma*beta**m*dx*xGrady(x)))
            m += 1

        if m<(converLevel-1):
            m_mem.append(m)
            x += beta**m*dx
            
            if verbosity != 'low':
                print('LS| x updated: '+str(x_mem[-1])+' -> '+str(x))
            x_mem.append(x)
        elif iloop == 0 and if_reversed_dx == False:
            print('LS| ***warning*** max Armijo order has been reached at the initial opt-step, searching direction is set\n'
                 +'                  reversed to find if possible to get lower value...')
            dx *= -1
            if_reversed_dx = True
            continue
        else:
            print('LS| ***warning*** max Armijo order has been reached, Armijo exit.')
            break
        iloop += 1
        if conservatism == False:
            iloop = 1
    
    return x

def armijo_step_revision(
    xfuncy, 
    xGrady, 
    x, 
    dx, 
    sigma = 0.2, 
    beta = 0.5, 
    m_thre = 20,
    verbosity = 'low'):
    # internal function, not be expected to be called directly from external function
    m = 0

    while (m <= m_thre) and (
        xfuncy(
         mop.matrix_plus(x, np.dot(beta**m, dx))
         ) 
         >= 
         (xfuncy(x) + np.dot(sigma*beta**m, np.dot(dx, xGrady(x))))
         ):
        
        if verbosity == 'high' or verbosity == 'debug':
            if verbosity == 'debug':
                print('OPT| Armijo: sigma = '+str(sigma)+', beta = '+str(beta))
            print('OPT| Armijo-asisted convergence method, step shrinking degree m = '+str(m))
        m += 1
    if m <= m_thre:
        return beta**m*dx
    else:
        print('OPT| ***warning*** Armijo method failed to find suitable stepsize, this may happen when\n'
             +'                   optimization task is already converged but inappropriate values are given\n'
             +'                   such as convergence threshold, Armijo beta, Armijo sigma.\n'
             +'     operation --> return stepsize = 0, to expire number of steps of the whole optimization.')
        return 0

def newton_1d(
    xfuncy, 
    x_ini, 
    dx_ini, 
    epsilon,
    conservatism = True,
    maxloop = 50,
    verbosity = 'low'):
    # note: epsilon here measures convergence threshold of variable value
    dx = dx_ini
    x_1 = x_ini
    x_2 = x_1 + dx
    x_3 = x_2 + dx
    # Newton assumes that it is already near minima, so curvature must be positive
    # make fitting quadratic function a*x**2 + b*x + c = f

    A = [
        [x_1**2, x_1, 1],
        [x_2**2, x_2, 1],
        [x_3**2, x_3, 1]
    ]
    B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
    if verbosity == 'debug':
        print('OPT| on-the-fly info: [type: initial info]\n'
             +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
             +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
             +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
             +'                      forward step: {}'.format(str(dx)))
    # solve for a, b, c
    para = np.linalg.solve(A, B)
    dx = -para[1]/2/para[0] - x_1
    iloop = 0
    while dx > float(epsilon) and iloop < maxloop:

        x_1 = -para[1]/2/para[0]
        x_2 = x_1 + dx
        x_3 = x_2 + dx
        A = [
            [x_1**2, x_1, 1],
            [x_2**2, x_2, 1],
            [x_3**2, x_3, 1]
        ]
        B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
        if verbosity == 'debug':
            print('OPT| on-the-fly info: [type: step info, step {}]\n'.format(str(iloop))
                +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
                +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
                +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
                +'                      forward step: {}'.format(str(dx)))
        para = np.linalg.solve(A, B)
        dx = -para[1]/2/para[0] - x_1
        iloop += 1

    return x_1

def steep_1d(
    xGrady, 
    x_ini, 
    acc_level, 
    shrink = True, 
    alpha = 0.1, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'):
    # original steep method
    # analytical expressions are needed
    x_mem = [x_ini, x_ini, x_ini]
    if shrink == False:
        alpha = 1
    dx = -alpha * xGrady(x_ini)

    iloop = 0
    while dx > 10**(-acc_level) and iloop < maxloop:
        x = x_mem[-1] + dx
        if verbosity == 'high' or verbosity == 'debug':
            print('OPT| steep method updates variable: '+str(x_mem[-1])+' -> '+str(x))
        x_mem = memoryUpdate(x_mem, x)
        dx = -alpha * xGrady(x)
        iloop += 1
        if conservatism == False:
            iloop = 1
    return x

# if function is multi-dimensional, its output is still 1-dimensional, but its derivative has the same 
# size as varaible, multi-dimensional. So for optimization method that not uses derivative, codes may 
# not need to be changed largely...

# from here, n-dimensional functions are supported

def armijo_steep(
    xfuncy,
    xGrady, 
    x_ini, 
    acc_level, 
    Armijo_sigma = 0.2, 
    Armijo_beta = 0.5, 
    max_Armijo = 20, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'
    ):
    x_mem = [x_ini, x_ini, x_ini]

    dx_0 = -xGrady(x_ini)
    dx = armijo_step_revision(
        xfuncy = xfuncy, 
        xGrady = xGrady,
        x = x_ini,
        dx = dx_0,
        sigma = Armijo_sigma,
        beta = Armijo_beta,
        m_thre = max_Armijo,
        verbosity = verbosity
        )

    iloop = 0
    while np.linalg.norm(dx) > 10**(-acc_level) and iloop < maxloop:
        x = mop.matrix_plus(x_mem[-1], dx)
        if verbosity == 'high' or verbosity == 'debug':
            print('OPT| steep method updates variable: '+str(x_mem[-1])+' -> '+str(x))
        x_mem = memoryUpdate(x_mem, x)
        dx_0 = -xGrady(x)
        dx = armijo_step_revision(
            xfuncy = xfuncy, 
            xGrady = xGrady,
            x = x_ini,
            dx = dx_0,
            sigma = Armijo_sigma,
            beta = Armijo_beta,
            m_thre = max_Armijo,
            verbosity = verbosity
            )
        iloop += 1
        if conservatism == False:
            iloop = 1
    return x

def newton(
    xfuncy, 
    x_ini, 
    dx_ini, 
    epsilon, 
    conservatism = True, 
    maxloop = 50, 
    verbosity = 'low'
    ):
    # multi-dimensional version of newton, also quadratic method
    
    ndim = len(x_ini)
    dx = dx_ini
    x_1 = x_ini
    x_2 = mop.matrix_plus(x_1, dx)
    x_3 = mop.matrix_plus(x_2, dx)
    x_1_sqr = mop.matrix_dot_power(x_1, 2)
    x_2_sqr = mop.matrix_dot_power(x_2, 2)
    x_3_sqr = mop.matrix_dot_power(x_3, 2)

    B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]

    if verbosity == 'debug':
        print('OPT| on-the-fly info: [type: initial info]\n'
             +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
             +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
             +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
             +'                      forward step: {}'.format(str(dx)))
    
    para = np.zeros((ndim, 3))

    for i in range(ndim):
        A = [
            [x_1_sqr[i], x_1[i], 1],
            [x_2_sqr[i], x_2[i], 1],
            [x_3_sqr[i], x_3[i], 1]
        ]
        para_this = np.linalg.solve(A, B)
        para[i][:] = para_this
        dx[i] = -para[i][1]/2/para[i][0] - x_1[i]

    iloop = 0

    while np.linalg.norm(dx) > float(epsilon) and iloop < maxloop:

        x_1 = np.dot(-0.5, mop.matrix_dot_division(para[:][1], para[:][0]))
        x_2 = mop.matrix_plus(x_1, dx)
        x_3 = mop.matrix_plus(x_2, dx)
        x_1_sqr = mop.matrix_dot_multiply(x_1, x_1)
        x_2_sqr = mop.matrix_dot_multiply(x_2, x_2)
        x_3_sqr = mop.matrix_dot_multiply(x_3, x_3)

        B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
        for i in range(ndim):
            A = [
                [x_1_sqr[i], x_1[i], 1],
                [x_2_sqr[i], x_2[i], 1],
                [x_3_sqr[i], x_3[i], 1]
            ]
            para_this = np.linalg.solve(A, B)
            para[i][:] = para_this
            dx[i] = -para[i][1]/2/para[i][0] - x_1[i]

        if verbosity == 'debug':
            print('OPT| on-the-fly info: [type: step info, step {}]\n'.format(str(iloop))
                +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
                +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
                +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
                +'                      forward step: {}'.format(str(dx)))

        iloop += 1

    return x_1

def armijo_steep_newton(
    xfuncy,
    xGrady, 
    x_ini, 
    dx_ini, 
    acc_level, 
    Armijo_sigma = 0.2, 
    Armijo_beta = 0.5, 
    max_Armijo = 20, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'
    ):
    epsilon = 10**(-acc_level)

    ndim = len(x_ini)
    dx = dx_ini
    x_1 = x_ini
    x_2 = mop.matrix_plus(x_1, dx)
    x_3 = mop.matrix_plus(x_2, dx)
    x_1_sqr = mop.matrix_dot_power(x_1, 2)
    x_2_sqr = mop.matrix_dot_power(x_2, 2)
    x_3_sqr = mop.matrix_dot_power(x_3, 2)

    B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]

    if verbosity == 'debug':
        print('OPT| on-the-fly info: [type: initial info]\n'
             +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
             +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
             +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
             +'                      forward step: {}'.format(str(dx)))
    
    para = np.zeros((ndim, 3))

    for i in range(ndim):
        A = [
            [x_1_sqr[i], x_1[i], 1],
            [x_2_sqr[i], x_2[i], 1],
            [x_3_sqr[i], x_3[i], 1]
        ]
        para_this = np.linalg.solve(A, B)
        para[i][:] = para_this
        dx[i] = -para[i][1]/2/para[i][0] - x_1[i]

    iloop = 0

    ifSteep = False
    steeplist = []
    while np.linalg.norm(dx) > float(epsilon) and iloop < maxloop:

        if ifSteep:
            x_1_prev = x_1
            x_1 = np.dot(-0.5, mop.matrix_dot_division(para[:][1], para[:][0]))
            for index in steeplist:
                x_1[index] = x_1_prev[index] + dx[index]
            ifSteep = False
            steeplist = []
        else:
            x_1 = np.dot(-0.5, mop.matrix_dot_division(para[:][1], para[:][0]))

        x_2 = mop.matrix_plus(x_1, dx)
        x_3 = mop.matrix_plus(x_2, dx)
        x_1_sqr = mop.matrix_dot_multiply(x_1, x_1)
        x_2_sqr = mop.matrix_dot_multiply(x_2, x_2)
        x_3_sqr = mop.matrix_dot_multiply(x_3, x_3)

        B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
        for i in range(ndim):
            A = [
                [x_1_sqr[i], x_1[i], 1],
                [x_2_sqr[i], x_2[i], 1],
                [x_3_sqr[i], x_3[i], 1]
            ]
            para_this = np.linalg.solve(A, B)
            para[i][:] = para_this
            # save parameters, anyway.
            if para[i][0] >= 0:
                dx[i] = -para[i][1]/2/para[i][0] - x_1[i]
            else:
                print('OPT| Armijo-steep-newton: negative curvature detected, optimization switches to Armijo-steep.')
                dx_0 = -xGrady(x_1[i])
                dx[i] = armijo_step_revision(
                        xfuncy = xfuncy, 
                        xGrady = xGrady,
                        x = x_ini,
                        dx = dx_0,
                        sigma = Armijo_sigma,
                        beta = Armijo_beta,
                        m_thre = max_Armijo,
                        verbosity = verbosity
                        )
                steeplist.append(i)
        if len(steeplist) > 0:
            ifSteep = True

        if verbosity == 'debug':
            print('OPT| on-the-fly info: [type: step info, step {}]\n'.format(str(iloop))
                +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
                +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
                +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
                +'                      forward step: {}'.format(str(dx)))

        iloop += 1

    return x_1
    