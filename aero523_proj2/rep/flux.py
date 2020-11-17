import numpy as np

def RoeFlux(Ul, Ur, n, return_test):
    gam = 1.4

    # Left-Side arguments
    rhol = Ul[0]; ul = Ul[1]/rhol; vl = Ul[2]/rhol
    pl = (gam - 1)*(Ul[3] - 0.5*rhol*(ul**2 + vl**2))
    Hl = (Ul[3] + pl)/rhol

    # Right-Side arguments
    rhor = Ur[0]; ur = Ur[1]/rhor; vr = Ur[2]/rhor; 
    pr = (gam - 1)*(Ur[3] - 0.5*rhor*(ur**2 + vr**2))
    Hr = (Ur[3] + pl)/rhor

    # Left and Right side fluxes
    FL = np.array([np.dot([rhol*ul, rhol*vl], n),np.dot([rhol*ul**2 + pl, rhol*ul*vl], n),np.dot([rhol*ul*vl, rhol*vl**2 + pl], n),np.dot([rhol*ul*Hl, rhol*vl*Hl], n)])
    FR = np.array([np.dot([rhor*ur, rhor*vr], n),np.dot([rhor*ur**2 + pr, rhor*ur*vr], n),np.dot([rhor*ur*vr, rhor*vr**2 + pr], n),np.dot([rhor*ur*Hr, rhor*vr*Hr], n)])
    
    RHS = ROE_Avg(ul,vl,ur,vr, rhol, rhor, Hl, Hr, pr, pl)
    F = 0.5*(FL + FR) - 0.5*RHS

    if return_test:
        return F, FL, FR
    else:
        return F

def ROE_Avg(ul,vl,ur,vr, rhol, rhor, Hl, Hr, Pr, Pl):
    vell = np.sqrt(ul**2 + vl**2); velr = np.sqrt(ur**2 + vr**2)
    url = (np.sqrt(rhol)*vell + np.sqrt(rhor)*velr)/(np.sqrt(rhol) + np.sqrt(rhor))
    hrl = (np.sqrt(rhol)*Hl + np.sqrt(rhor)*Hr)/(np.sqrt(rhol) + np.sqrt(rhor))
    rhorl = np.sqrt(rhor*rhol)
    arl = np.sqrt((1.4 - 1.0)*(hrl - 0.5*url**2))

    ls = np.array([url, url+arl, url-arl])
    ws = np.array([rhor-rhol - (Pr-Pl)/arl**2, ur-ul + (Pr-Pl)/(rhorl*arl), ur-ul - (Pr-Pl)/(rhorl*arl)])
    es = np.array([[1, 1, 1], [1, 1, 1], [url, url + arl, url-arl], [1/2*url**2, hrl+url*arl, hrl - url*arl]])
    es[:,1] *= rhorl/(2*arl); es[:,2] *= -rhorl/(2*arl)

    RHS = 0
    for i in range(3):
        if abs(ls[i]) < 0.1*arl:
            phil = (ls[i]**2 + (0.1*arl)**2)/(2*0.1*arl)
        else:
            phil = abs(ls[i])
        RHS += 0.5*phil*es[:,i]*ws[i]
    
    return RHS
