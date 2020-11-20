import numpy as np
from numpy import linalg as LA

def RoeFlux(Ul, Ur, n, return_test):
    gam = 1.4

    # Left side arguments
    rhol = Ul[0]; ul = Ul[1]/rhol; vl = Ul[2]/rhol; rhoEl = Ul[3]
    pl = (gam-1)*(rhoEl-0.5*rhol*(ul**2 + vl**2))
    Hl = (rhoEl + pl)/rhol

    # Right side arguments
    rhor = Ur[0]; ur = Ur[1]/rhor; vr = Ur[2]/rhor; rhoEr = Ur[3]
    pr = (gam-1)*(rhoEr-0.5*rhor*(ur**2 + vr**2))
    Hr = (rhoEr + pr)/rhor

    # Left and Right side fluxes
    FL = np.array([np.dot([Ul[1],Ul[2]], n), np.dot([Ul[1]*ul+pl, Ul[2]*ul],n), np.dot([Ul[1]*vl, Ul[2]*vl+pl],n), Hl*np.dot([Ul[1],Ul[2]],n)])
    FR = np.array([np.dot([Ur[1],Ur[2]], n), np.dot([Ur[1]*ur+pr, Ur[2]*ur],n), np.dot([Ur[1]*vr, Ur[2]*vr+pr],n), Hr*np.dot([Ur[1],Ur[2]],n)])
    
    # Roe-Averages
    RHS, ls = ROE_Avg(ul,vl,rhol,Hl,rhoEl, ur,vr,rhor,Hr,rhoEr, n)
    F = 0.5*(FL + FR) - 0.5*RHS
   
    if return_test:
        return F, FL, FR
    else:
        return F, ls

def ROE_Avg(ul,vl,rhol,Hl,rhoEl, ur,vr,rhor,Hr,rhoEr, n):
    gam = 1.4
    vell = np.array([ul, vl]); velr = np.array([ur, vr])

    # Calculating Roe average
    v = (np.sqrt(rhol)*vell + np.sqrt(rhor)*velr)/(np.sqrt(rhol) + np.sqrt(rhor))
    H = (np.sqrt(rhol)*Hl + np.sqrt(rhor)*Hr)/(np.sqrt(rhol) + np.sqrt(rhor))

    # Calculating eigenvalues
    q = LA.norm(v)
    c = np.sqrt((gam-1.0)*(H - 0.5*q**2))
    u = np.dot(v, n)
    ls = np.array([u+c, u-c, u])

    # Apply the entropy fix
    for i in range(3):
        if abs(ls[i]) < 0.1*c:
            ls[i] = ((0.1*c)**2 + ls[i]**2)/(2*0.1*c)
    ls = abs(ls)

    delrho = rhor - rhol; delmo = np.array([rhor*ur - rhol*ul, rhor*vr - rhol*vl]); dele = rhoEr - rhoEl
    s1 = 0.5*(abs(ls[0]) + abs(ls[1])); s2 = 0.5*(abs(ls[0]) - abs(ls[1]))
    G1 = (gam-1.0)*(0.5*q**2*delrho - np.dot(v, delmo) + dele); G2 = -u*delrho + np.dot(delmo, n)
    C1 = G1*(c**-2)*(s1 - abs(ls[2])) + G2*(c**-1)*s2; C2 = G1*(c**-1)*s2 + (s1 - abs(ls[2]))*G2

    RHS = np.array([ls[2]*delrho+C1, ls[2]*delmo[0]+C1*v[0]+C2*n[0], ls[2]*delmo[1]+C1*v[1]+C2*n[1], ls[2]*dele+C1*H+C2*u]) 

    return RHS, max(ls)

def Rus(Ul,pl, Ur,pr):
    gam = 1.4
    cl = np.sqrt(gam*pl/Ul[0])
    cr = np.sqrt(gam*pr/Ur[0])

    s = np.array([abs(Ul[2]/Ul[0]) + cl, abs(Ul[1]/Ul[0]) + cl, abs(Ur[2]/Ur[0] + cr), abs(Ur[1]/Ur[0] + cr)])
    
    return max(s)

