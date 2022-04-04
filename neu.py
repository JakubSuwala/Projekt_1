# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:40:48 2022

@author: USER
"""

def s_A_z2neu(s, A, z):
    #A = np.deg2rad(A)
    #z = np.deg2rad(z)
    n = s*np.sin(z)*np.cos(A)
    e = s*np.sin(z)*np.sin(A)
    u = s*np.cos(z)
    return n, e, u