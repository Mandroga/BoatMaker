import numpy as np
from scipy.optimize import minimize

class BaseMastroSolver:
    def __init__(self, rm, ym, u0, LS, R = 100, er = 10**-3):
        self.rm = rm
        self.ym = ym
        self.u0 = u0
        b3 = ym-rm-u0

        index = np.argmin(abs(LS.y4-b3))

        y1 = LS.y4[index:]
        x1R = LS.x4[index:]
        x1L = -x1R

        t1 = np.linspace(0,2*np.pi,R)
        x2 = rm*np.cos(t1)
        y2 = rm*np.sin(t1)+ym

        y3 = np.array([ym-rm-u0]*2)
        x3 = np.linspace(x1L[0], x1R[0], 2)

        self.A = np.trapz(x1R,y1)
        self.Cy = np.trapz(x1R * y1, y1) / self.A
        dr = 0.004
        d_m = 700
        self.P = self.A * dr * d_m

        self.lm = x1R[0]-x1L[0]
        self.b3 = b3
        self.plot = [[x1R,x1L,x2,x3],[y1,y1,y2,y3]]


