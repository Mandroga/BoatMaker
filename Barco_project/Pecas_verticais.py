import numpy as np

class VerticaisSolver:
    def __init__(self, oi, l2, h2, ym, rm, u0, LS, BMS):
        l4 = LS.l4
        o = LS.op(ym-rm-u0)
        b3 = BMS.b3
        lm4 = BMS.lm
        index = np.argmin(abs(LS.y3-b3))
        lm3 = 2*LS.x3[index]
        c1 = l4/2-l2/2
        c2 = lm4/2-lm3/2


        x1 = np.linspace(-l4/2, l4/2,2)
        x2 = np.linspace(-l2/2, l2/2,2)
        x3 = np.linspace(l2/2,l4/2,2)

        x4 = np.linspace(-lm4/2, lm4/2,2)
        x5 = np.linspace(-lm3/2,lm3/2,2)
        x6 = np.linspace(lm3/2,lm4/2,2)

        def side(x):
            return (h2/c1)*(x-l2/2)
        def side2(x):
            return (h2/c2)*(x-lm3/2)+h2

        self.A = (h2/2)*np.array([l4+l2,lm4+lm3])
        dr = 0.004
        d_m = 700
        self.P = self.A*dr*d_m
        self.Cy = np.array([LS.b2,BMS.b3])
        self.Cz = np.array([])

        self.plot = [[x1,x2,x3,-x3, x4,x5,x6,-x6],[[h2]*2,[0]*2,side(x3),side(x3), [h2*2]*2,[h2]*2, side2(x6), side2(x6)]]

