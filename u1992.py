# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:44:50 2022

@author: USER
"""

def u1992(fi,lam,a,e2):
    m_0 = 0.9993
    N = a/(m.sqrt(1-e2 * np.sin(fi)**2))
    t = np.tan(fi)
    n2 = e_2 * np.cos(lam)**2
    lam_0 = m.radians(19)
    l = lam - lam_0
    
    A_0 = 1 - (e2/4) - (3*(e2**2))/64 - (5*(e2**3))/256
    A_2 = 3/8 * (e2 + ((e2**2)/4) + ((15*e2**3)/128))
    A_4 = 15/256 * (e2**2 + (3*(e2**3))/4)
    A_6 = (35*(e2**3))/3072
    
    sigma = a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
    
    x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
    y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

    x92 = round(x * m_0 - 5300000, 3)
    y92 = round(y * m_0 + 500000, 3)   
    
    return x92, y92 