import numpy as np

def Fluxes(ul, ur, flux_type):
    gam = 1.4
    
    Pl = (gam-1)*(ul[2] - 0.5*ul[1]**2/ul[0])
    Fl = np.array([ul[1], ul[1]**2/ul[0] + Pl, ul[1]*(ul[2]/ul[0] + Pl/ul[0])])
    Pr = (gam-1)*(ur[2] - 0.5*ur[1]**2/ur[0])
    Fr = np.array([ur[1], ur[1]**2/ur[0] + Pr, ur[1]*(ur[2]/ur[0] + Pr/ur[0])])

    cl = np.sqrt(gam*Pl/ul[0])
    cr = np.sqrt(gam*Pr/ur[0])

    if flux_type == 'Rus':
        s = np.array([abs(ul[1]/ul[0]) + cl, abs(ur[1]/ur[0] + cr)])
        
        flux = 0.5*(Fl + Fr) - 0.5*max(s)*(ur - ul)
    if flux_type == 'HLLE':
        smax = max(np.array([abs(ul[1]/ul[0]) + cl, abs(ur[1]/ur[0] + cr)]))
        smin = min(np.array([abs(ul[1]/ul[0]) - cl, abs(ur[1]/ur[0] - cr)]))
        
        flux = 0.5*(Fl + Fr) - 0.5*((smax + smin)/(smax - smin))*(Fr - Fl) + ((smax * smin)/(smax - smin))*(ur - ul)
    
    return np.transpose(flux)

def getIC(ul, ur, x):
    u0 = np.zeros((3, x.shape[0]))
    for i in range(x.shape[0]):
        if x[i] < 0.5:
            u0[:,i] = ul
        else:
            u0[:,i] = ur
    
    return u0

def l2err(R):
    err = np.zeros((3,1))
    for i in range(R.shape[0]):
        err[:,0] += (0.0 - R[:,i])**2
    err = np.sqrt(1/R.shape[0] * err)
    
    return np.transpose(err)

def solve(ul0, ur0, N, T, verify):
    L = 1; gam = 1.4
    x = np.linspace(0, L, num=N+1, endpoint=True)
    xc = 0.5*(x[0:N] + x[1:N+1])
    dx = xc[1] - xc[0]
    
    Pl = (gam-1)*(ul0[2] - 0.5*ul0[1]**2/ul0[0])
    Pr = (gam-1)*(ur0[2] - 0.5*ur0[1]**2/ur0[0])
    cl = np.sqrt(gam*Pl/ul0[0])
    cr = np.sqrt(gam*Pr/ur0[0])
    a = max(np.array([cl, cr]))
    dt = 0.5*dx/a; Nt = int(np.ceil(T/dt)); dt = T/Nt
    #Nt = 1000; dt = T/Nt # Un-comment if verifying
    
    u0 = getIC(ul0, ur0, xc)
    uRus = u0.copy(); RRus = uRus.copy()
    uHLLE = u0.copy(); RHLLE = uHLLE.copy()
    Rus_resid = np.zeros((3,Nt)); HLLE_resid = Rus_resid.copy()
    for n in range(Nt):
        RRus *= 0; RHLLE *= 0
        for j in range(N+1):
            ul = uRus[:,j-1] if (j > 0) else ul0
            ur = uRus[:,j  ] if (j < N) else ur0
            FRus = Fluxes(ul,ur, 'Rus')

            ul = uHLLE[:,j-1] if (j > 0) else ul0
            ur = uHLLE[:,j  ] if (j < N) else ur0
            FHLLE = Fluxes(ul,ur, 'HLLE')
            if (j > 0): 
                RRus[:,j-1] += FRus
                RHLLE[:,j-1] += FHLLE
            if (j < N):
                RRus[:,j  ] -= FRus
                RHLLE[:,j] -= FHLLE

        if verify:
            Rus_resid[:,n] = l2err(RRus)
            HLLE_resid[:,n] = l2err(RHLLE)

        uRus -= dt/dx * RRus
        uHLLE -= dt/dx * RHLLE
        
    return xc, uRus, uHLLE, Rus_resid, HLLE_resid
    
    