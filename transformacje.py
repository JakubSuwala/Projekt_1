from math import sin, cos, sqrt, atan, atan2, degrees, radians
import math
import numpy as np

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przyblizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            
            
    def filamh2XYZ(self, fi, lam, h):
        '''
        funkcja dokonuje transformacji współrzędnych geodezyjnych na ortokartezjańskie.
        '''
        fi = math.radians(fi)
        lam = math.radians(lam)
        N = self.a / math.sqrt(1 - self.ecc2 * math.sin(fi)**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N * (1 - self.ecc2) + h) * math.sin(fi)
        return(X, Y, Z)



    def neu(self, Fia, La, Fib, Lb, h):
        '''
        funkcja dokonuje transformacji współrzędnych do układu neu.
        '''
    
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
    
    
    def u2000(self, fi, lam):
        '''
        funkcja dokonuje transformacji współrzędnych do układu 2000.
        '''
        e_2 = self.ecc2 / (1 - self.ecc2)
        m_0 = 0.999923
        N = self.a / (math.sqrt(1 - self.ecc2 * np.sin(fi) ** 2))
        t = np.tan(fi)
        n2 = e_2 * np.cos(lam) ** 2
        # lam = math.degrees(lam)
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
        lam = math.radians(lam)
        lam_0 = math.radians(lam_0)
        l = lam - lam_0
        A_0 = 1 - (self.ecc2 / 4) - (3 * (self.ecc2 ** 2)) / 64 - (5 * (self.ecc2 ** 3)) / 256
        A_2 = 3 / 8 * (self.ecc2 + ((self.ecc2 ** 2) / 4) + ((15 * self.ecc2 ** 3) / 128))
        A_4 = 15 / 256 * (self.ecc2 ** 2 + (3 * (self.ecc2 ** 3)) / 4)
        A_6 = (35 * (self.ecc2 ** 3)) / 3072
        sigma = self.a * ((A_0 * fi) - (A_2 * np.sin(2 * fi)) + (A_4 * np.sin(4 * fi)) - (A_6 * np.sin(6 * fi)))
        x = sigma + ((l ** 2) / 2) * (N * np.sin(fi) * np.cos(fi)) * (
                    1 + ((l ** 2) / 12) * ((np.cos(fi)) ** 2) * (5 - t ** 2 + 9 * n2 + (4 * n2 ** 2)) + ((l ** 4) / 360) * (
                        (np.cos(fi)) ** 4) * (61 - (58 * (t ** 2)) + (t ** 4) + (270 * n2) - (330 * n2 * (t ** 2))))
        y = l * (N * np.cos(fi)) * (1 + ((((l ** 2) / 6) * (np.cos(fi)) ** 2) * (1 - (t ** 2) + n2)) + (
                    ((l ** 4) / (120)) * (np.cos(fi) ** 4)) * (
                                                5 - (18 * (t ** 2)) + (t ** 4) + (14 * n2) - (58 * n2 * (t ** 2))))
        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s * 1000000) + 500000, 3)
    
        return x00, y00


    def u1992(self, fi, lam):
        '''
        funkcja dokonuje transformacji współrzędnych do układu 1992.
        '''
        m_0 = 0.9993
        e_2 = self.ecc2 / (1 - self.ecc2)
        N = self.a/(math.sqrt(1-self.ecc2 * np.sin(fi)**2))
        t = np.tan(fi)
        n2 = e_2 * np.cos(lam)**2
        lam_0 = math.radians(19)
        l = lam - lam_0
        
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        return x92, y92 
    
    
    def odl_2D_3D (self, X, Y, Z):
        '''
        funkcja dokonuje obliczeń odległosci 2D oraz 3D dla wartosci współrzędnych ortokartezjańskich.
        '''
        odl_2D = math.sqrt(X**2 + Y**2)
        odl_3D = math.sqrt(X**2 + Y**2 + Z**2)
        
        return odl_2D, odl_3D

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    
    X, Y, Z = geo.filamh2XYZ(phi, lam, h)
    # print(X, Y, Z)
        
    n, e, u = geo.neu(phi, lam, 50, 20, h) # przykladowe dane do testu
    # print(n, e, u)
    
    x00, y00 = geo.u2000(phi, lam)
    # print(x00, y00)
    
    x92, y92 = geo.u1992(phi, lam)
    # print(x92, y92)
    
    odl_2D, odl_3D = geo.odl_2D_3D(X, Y, Z)
    # print(odl_2D, odl_3D)
    
    
    
    