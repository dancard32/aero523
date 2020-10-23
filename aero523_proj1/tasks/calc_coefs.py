def calc_cp(ylin, phi):
    import numpy as np

    h = ylin[0,0] - ylin[1,0]
    num = max(np.shape(ylin)) - 1
    cp = np.zeros(num + 1)
    
    for i in range(num+1):
        u = (-3/2*phi[num,i] + 2*phi[num-1,i] - 1/2*phi[num-2,i])*h**-1
        cp[i] = 1 - u**2

    return cp

def calc_cl(ylin, cp):
    import numpy as np

    h = ylin[0,0] - ylin[1,0]
    num = max(np.shape(ylin)) - 1
    cl = 0
    
    for i in range(1, num):
        cl += h*cp[i]
    cl += h/2*(cp[0] + cp[num])
    cl *= 1/2

    return cl