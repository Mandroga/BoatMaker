import numpy as np
import matplotlib.pyplot as plt
from Base_Barco_backup import BaseSolver
from Tools import ComputationalMath as CM
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize

def Lateral_Solver(BS, h2, oi ,of, er=10**-5, R=1000):

    r3 = BS[0]
    a = BS[1]
    i1 = BS[2]
    i2 = BS[3]
    b1 = BS[4]
    b2 = BS[5]

    phi = 2*np.pi*(1 - (np.tan(oi)/np.sqrt(1+np.tan(oi)*np.tan(oi))))
    
    r4max = (h2 * np.tan(of) + r3)
    b0 = np.sqrt(r4max * r4max - a * a)

    def op(y):
        return ((of-oi)/(b0-b2)) * (y-b2) + oi
    
    ycritmax = (np.pi/2-oi-0.01)*( (b0-b2)/(of-oi) )+b2
    ycritmin = (-oi)*( (b0-b2)/(of-oi) )+b2
    
    def r4(y):
        return h2*np.tan(op(y))+r3
    
    r1 = r3/(1-phi/(2*np.pi))
    
    def r2(y):
        return r4(y)/(1-phi/(2*np.pi))
    
    l3 = 2*(-a+r4(0))
    l4 = 2*(-a+np.sqrt(r4(b2)*r4(b2)-b2*b2))
    
    i3 = np.arctan(b2/(a+l4/2))
    i4 = np.arctan(b0/a)
    
    ii1 = (r3/r1)*i1
    ii2 = (r3/r1)*i2
    
    def r1_line(y):
        return -(a - np.sqrt(r1*r1-y*y))
    
    def r2_line(y):
        return -(a-np.sqrt(r2(y)*r2(y)-y*y))
    
    def r3_line(y):
        return -(a - np.sqrt(r3*r3-y*y))
    
    def r4_line(y):
        return  -(a - np.sqrt(r4(y)*r4(y)-y*y))
    
    k1 = h2/np.cos(oi)
    k2 = h2/np.cos(of)
    
    
    def y2isolver(y):
        dx = r1*np.cos(ii1)-a - r2_line(y)
        dy = r1*np.sin(ii1) - y
        return abs(k1 - np.sqrt(dx*dx + dy*dy))
    
    y1i = r1*np.sin(ii1)
    y1f = r1*np.sin(ii2)
    
    y2i = minimize(y2isolver, 0).x[0]
    y3 = np.linspace(b2,b1,R)
    y4 = np.linspace(b2,b0,R)
    y1 = np.linspace(r1*np.sin(ii1),r1*np.sin(ii2),R)
    
    L4 = CM.Lenght(r4_line(y4),y4)
    #L2 = CM.Lenght(r2_line(y2), y2)[0]

    def y2fsolver():
        y2f = y1f
        y = np.linspace(y2i, y2f,100)
        L = CM.Lenght(r2_line(y),y)
    
        while L<L4:
            p0 = [r2_line(y2f), y2f]
            y2f += er
            p1 = [r2_line(y2f), y2f]
            L += CM.Lenght([p0[0], p1[0]], [p0[1], p1[1]])
        return y2f
    
    y2f = y2fsolver()
    y2 = np.linspace(y2i, y2f, R)

    x5 = np.linspace(r1_line(y1i),r2_line(y2i),2)
    x6 = np.linspace(r1_line(y1f),r2_line(y2f),2)
    y5 = np.linspace(y1i, y2i, 2)
    y6 = np.linspace(y1f, y2f, 2)
    return [[r1_line(y1),r2_line(y2),r3_line(y3),r4_line(y4), x5, x6], [y1,y2,y3,y4,y5,y6]]