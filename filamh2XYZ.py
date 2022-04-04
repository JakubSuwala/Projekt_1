# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:38:25 2022

@author: USER
"""

def filamh2XYZ(fi, lam, h, a, e2):
    N = a/math.sqrt(1-e2*math.sin(fi)**2)
    X = (N + h) * math.cos(fi) * math.cos(lam)
    Y = (N + h) * math.cos(fi) * math.sin(lam)
    Z = (N*(1-e2) + h) * math.sin(fi)
    return(X, Y, Z)