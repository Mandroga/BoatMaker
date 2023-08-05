import numpy as np
import matplotlib.pyplot as plt


class BaseSolver:
    def __init__(self, l1, l2, h1, R = 100):
        b = (l2*l2)/(8*h1)-(l1*l2)/(4*h1)+h1/2
        c = l2/(2*h1)
        d = 2*c*b - l1
        e = b*b + (l1*l1)/4
        f = c*c

        r = (-d-np.sqrt(d*d-4*f*e))/(2*f)
        a = r - l1/2

        b1 = np.sqrt(r*r-a*a)
        b2 = b1 - h1

        i2 = np.arctan(b2/(a+l2/2))
        i1 = np.arctan(b1/a)



        def Boat_BaseL(y):
            return a - np.sqrt(r * r - y * y)

        def Boat_BaseR(y):
            return - a + np.sqrt(r * r - y * y)

        y = np.linspace(b2, b1, R)
        x = Boat_BaseR(y)

        self.A = np.trapz(x,y-b2)
        dr = 0.004
        d_m = 700
        self.P = self.A * dr * d_m
        self.Cy = np.trapz(x*y,y)/self.A

        self.r = r
        self.a = a
        self.i1 = i1
        self.i2 = i2
        self.b1 = b1
        self.b2 = b2
        self.l1 = l1
        self.l2 = l2
        self.plot = [[Boat_BaseL(y),Boat_BaseR(y),np.linspace(-l2/2,l2/2,3)],[y,y,[b2]*3]]

