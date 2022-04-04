# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:49:19 2022

@author: USER
"""

def odl_2D_3D (X, Y, Z):
    odl_2D = math.sqrt(X**2 + Y**2)
    odl_3D = math.sqrt(X**2 + Y**2 + Z**2)
    return odl_2D, odl_3D