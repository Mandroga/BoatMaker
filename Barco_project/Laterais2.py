import numpy as np
import matplotlib.pyplot as plt
from Base import BaseSolver
from Tools import ComputationalMath as CM
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import curve_fit

class LateralSolver:
    def __init__(self, BS, h2, oi, of, er=10 ** -5, R=1000):
            r3 = BS.r
            a = BS.a
            i1 = BS.i1
            i2 = BS.i2
            b1 = BS.b1
            b2 = BS.b2
            l2 = BS.l2

            v1 = h2 * np.tan(of)
            b0 = b1 + v1

            i3 = np.arctan(b0 / a)

            c2 = np.sqrt(a * a + b0 * b0) - r3
            ofp = np.arctan(c2 / h2)

            k2 = h2 / np.cos(of)
            k1 = h2 / np.cos(oi)

            def op(i):
                    return ((ofp - oi) / (i3 - i2)) * (i - i2) + oi

            def r4(i):
                    return r3 + h2 * np.tan(op(i))

            i4 = minimize(lambda i: abs(r4(i)*np.sin(i)-b2),i2).x[0]

            t3 = np.linspace(i2, i1, R)
            x3 = r3 * np.cos(t3) - a
            y3 = r3 * np.sin(t3)
            L3 = CM.Lenght(x3, y3)

            t4 = np.linspace(i4, i3, R)
            x4 = r4(t4) * np.cos(t4) - a
            y4 = r4(t4) * np.sin(t4)
            L4 = CM.Lenght(x4, y4)

            self.l4 = 2*x4[0]

            r1 = r3 / np.sin(oi)

            ii1 = np.sin(oi) * i1
            ii2 = np.sin(oi) * i2
            ii3 = np.sin(oi) * i3
            ii4 = np.sin(oi) * i4

            t1 = np.linspace(ii2, ii1, R)
            x1 = r1 * np.cos(t1) - r1 + 1
            y1 = r1 * np.sin(t1)
            L1 = CM.Lenght(x1, y1)

            di = 0

            def kiSolver(i):
                kt = np.linspace(i4, i, R)
                k = h2 / np.cos(op(kt))

                t2 = np.linspace(ii4, i * np.sin(oi), R)
                r2 = r1 + k
                x2 = r2 * np.cos(t2) - r1 + 1
                y2 = r2 * np.sin(t2)
                L2 = CM.Lenght(x2, y2)
                k2m = CM.Lenght([x1[-1], x2[-1]], [y1[-1], y2[-1]])
                return abs(L4-L2)+abs(k2-k2m)

            ki = minimize(kiSolver,i3).x[0]
            kii = ki * np.sin(oi)

            print(f"ki = {ki}, kii = {kii}")

            kt = np.linspace(i4, ki, R)
            k = h2 / np.cos(op(kt))
            self.k = k
            self.kt = kt

            t2 = np.linspace(ii4, kii, R)
            self.t2 = t2
            r2 = r1 + k
            x2 = r2 * np.cos(t2) - r1 + 1
            y2 = r2 * np.sin(t2)
            L2 = CM.Lenght(x2, y2)

            print(f"L4 = {L4}, L2 = {L2}")

            k2m = CM.Lenght([x1[-1], x2[-1]], [y1[-1], y2[-1]])
            k1m = CM.Lenght([x1[0], x2[0]], [y1[0],  y2[0]])

            x5 = np.array([x1[0],x2[0]])
            y5 = np.array([y1[0],y2[0]])


            #print(f"L2-L4 = {L2 - L4}, L1-L3 = {L1 - L3}")
            #print(f"k2 - k2m = {k2 - k2m}")

            def kfit(i,a,b,c):
                return a/(i-b)+c

            fit_params, covariance = curve_fit(kfit,t2,k,[1,-5,0.13])
            self.fit_params = fit_params

            def r2(i):
                a, b, c = fit_params
                return r1 + kfit(i,a,b,c)

            def f(i):
                _x1 = x1[-1]+r1-1
                _x2 = x2[-1]+r1-1
                a = (y2[-1]-y1[-1])/(_x2-_x1)
                b = _x1
                c = y1[-1]
                return (c-a*b)/(np.sin(i)-a*np.cos(i))

            #Area of a single piece
            A = np.trapz(0.5*(r2(t1)*r2(t1)-r1*r1),t1)

            tA = np.linspace(ii1, kii, R)
            A += np.trapz(0.5*(r2(tA)*r2(tA)-f(tA)*f(tA)),tA)

            x6 = f(tA)*np.cos(tA)-r1+1
            y6 = f(tA)*np.sin(tA)

            plt.show()
            self.r3 = r3
            self.l2 = l2
            self.a = a
            self.h2 = h2
            self.b1 = b1
            self.oi = oi
            self.of = of
            self.ofp = ofp
            self.i1 = i1
            self.i2 = i2
            self.i3 = i3
            self.i4 = i4
            self.ki = ki
            self.kii = kii
            self.ii4 = ii4
            self.b0 = b0
            self.b2 = b2
            self.y4 = y4
            self.x4 = x4
            self.y3 = y3
            self.x3 = x3
            self.A = A
            print(f"A = {A}")
            self.plot = [[x1, x2, x5, x6, x3, x4], [y1, y2, y5, y6, y3, y4]]

    def op(self,i):
            return ((self.ofp - self.oi) / (self.i3 - self.i2)) * (i - self.i2) + self.oi

    def rh(self,i, hy):
        return self.r3 + hy * np.tan(self.op(i))

    def Impulsion(self, hyf, yp, R=100):
        hy = 0
        dh = hyf/R
        i_i = self.i2
        A = [0]*(len(yp)+1)
        while hy < hyf:
            i_i = minimize(lambda i: abs(self.rh(i,hy)*np.sin(i)-self.b2), i_i).x[0]
            i_f = np.arctan( (self.b1+hy*np.tan(self.of))/self.a)
            t = np.linspace(i_i, i_f, R)
            x = self.rh(t,hy) * np.cos(t) - self.a
            y = self.rh(t,hy) * np.sin(t)

            II = [0]
            for val in yp:
                II += [np.argmin(abs(y-val))]
            II += [-1]
            #II vai ter len 4
            #yp len 2
            ii = 0
            while ii < len(yp)+1:
                A[ii] += 2*np.trapz(x[II[ii]:II[ii+1]],y[II[ii]:II[ii+1]])
                ii += 1
            hy += dh
        A = np.array(A)
        V = A * dh
        #print(f" Volume = {V}")
        I = 997*V
        print(f"Impulsion in Kg = {I}")

    def r4(self,i):
            return self.r3 + self.h2 * np.tan(self.op(i))

    def kfit(self, i, a, b, c):
            return a / (i - b) + c