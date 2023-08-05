import numpy as np
import matplotlib.pyplot as plt


def BaseSolver(l1, l2, h1):
    b = (l2*l2)/(8*h1)-(l1*l2)/(4*h1)+h1/2
    c = l2/(2*h1)
    d = 2*c*b - l1
    e = b*b + (l1*l1)/4
    f = c*c

    r = (-d-np.sqrt(d*d-4*f*e))/(2*f)
    a = r - l1/2

    b1 = np.sqrt(r*r-a*a)
    b2 = b1 - h1

    i1 = np.arctan(b2/(a+l2/2))
    i2 = np.arctan(b1/a)

    y = np.linspace(b2, b1, 20)
    x = np.linspace(-l1, l1, 20)

    def Boat_BaseL(y):
        return a - np.sqrt(r * r - y * y)

    def Boat_BaseR(y):
        return - a + np.sqrt(r * r - y * y)

    A = 2*np.trapz(y-b2,Boat_BaseR(y-b2))

    return [[r, a, i1, i2, b1, b2, A],[[Boat_BaseL(y),Boat_BaseR(y),np.linspace(-l2/2,l2/2,3)],[y,y,[b2]*3]]]

