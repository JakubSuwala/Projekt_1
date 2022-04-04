# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:40:48 2022

@author: USER
"""
    def neu(self, Fia, La, Fib, Lb, h):
    
        Fia = radians(Fia)
        Fib = radians(Fib)
        La = radians(La)
        Lb = radians(Lb)

        dL = Lb - La
        Ua = math.atan((1-self.flat) * math.tan(Fia))
        Ub = math.atan((1-self.flat) * math.tan(Fib))
        Ln = dL
        while True:
            L = Ln
            sinS = math.sqrt((math.cos(Ub) * math.sin(L))**2 + (math.cos(Ua) * math.sin(Ub) - math.sin(Ua) * math.cos(Ub) * math.cos(L))**2)
            cosS = (math.sin(Ua) * math.sin(Ub) + math.cos(Ua) * math.cos(Ub) * math.cos(L))
            Si = math.atan2(sinS,cosS)
            sinA = (math.cos(Ua) * math.cos(Ub) * math.sin(L)) / sinS
            cos2A = 1 - sinA**2
            cos2S = cosS - (2*math.sin(Ua)*math.sin(Ub))/cos2A
            C = (self.flat/16) * cos2A * (4 + self.flat * (4 - 3*cos2A))
            Ln = dL + (1 - C) * self.flat * sinA*(Si + C * sinS * (cos2S + C*cosS*(-1+2*(cos2S**2))))
            if abs(Ln-L) < (math.radians(0.00001/3600)):
                L = Ln
                break
        u2 = ((self.a**2 - self.b**2)/(self.b**2)) * cos2A
        A = 1 + (u2/16384) * (4096 + u2*(-768 + u2*(320 - 175*u2)))
        B = (u2/1024) * (256 + u2*(-128 + u2*(74-47*u2)))
        dSi = B * sinS * (cos2S + (0.25) * B * (cosS*(-1+2*(cos2S**2)) -(B/6)*cos2S*(-3 + 4*(sinS**2))*(-3 + 4*(cos2S**2))))
        Sab = self.b * A * (Si - dSi)
        Aab = math.atan2(math.cos(Ub) * math.sin(L), math.cos(Ua) * math.sin(Ub) - math.sin(Ua) * math.cos(Ub) * math.cos(L))
        Aba = math.atan2(math.cos(Ua) * math.sin(L), -math.sin(Ua) * math.cos(Ub) + math.cos(Ua) * math.sin(Ub) * math.cos(L)) + math.pi
        if Aab<0: Aab = Aab + 2*math.pi
        if Aba<0: Aab = Aba + 2*math.pi
    
        z = 90*math.pi/180 - h
        s = Sab
        A = Aab
    
        #A = np.deg2rad(A)
        #z = np.deg2rad(z)
        n = s*np.sin(z)*np.cos(A)
        e = s*np.sin(z)*np.sin(A)
        u = s*np.cos(z)
    
        return n, e, u