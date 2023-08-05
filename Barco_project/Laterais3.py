import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import json
from Base import BaseSolver

class LateralSolver:
    def __init__(self):
        print("Class Constructed")
    def Build(self, BS, h2, oi, of, er=10 ** -5, R=100):
        r3 = BS.r
        self.r3 = r3
        self.a = BS.a
        self.h2 = h2
        self.R = R

        v1 = h2 * np.tan(of)
        b0 = BS.b1 + v1

        self.b0 = b0
        self.b1 = BS.b1
        self.b2 = BS.b2

        c2 = np.sqrt(BS.a * BS.a + b0 * b0) - r3
        ofp = np.arctan(c2 / h2)

        self.oi = oi
        self.of = of
        self.ofp = ofp

        i3 = np.arctan(b0 / BS.a)

        self.i1 = BS.i1
        self.i2 = BS.i2
        self.i3 = i3
        self.i4 = minimize(lambda i: abs(self.rh(i, h2) * np.sin(i) - BS.b2), BS.i2).x[0]

        self.l4 = 2 * self.rh(self.i4,h2) * np.cos(self.i4) - BS.a

        self.r1 = r3 / np.sin(oi)

        self.ii1 = np.sin(oi) * BS.i1
        self.ii2 = np.sin(oi) * BS.i2
        self.ii4 = np.sin(oi) * self.i4

        self.ki, self.kii = self.kiSolver()

    def L(self,y,x):
        return np.trapz(y,x)

    def op(self,i):
            return ((self.ofp - self.oi) / (self.i3 - self.i2)) * (i - self.i2) + self.oi

    def rh(self,i, hy):
        return self.r3 + hy * np.tan(self.op(i))


    def r1_line(self):
        t = np.linspace(self.ii2, self.ii1, self.R)
        x = self.r1 * np.cos(t)
        y = self.r1 * np.sin(t)
        return x,y

    def r2(self, i):
        return self.r1 + self.h2 / np.cos(self.op(i))

    def r2_line(self):
        kt = np.linspace(self.i4, self.ki, self.R)

        t = np.linspace(self.ii4, self.kii, self.R)
        r2 = self.r2(kt)
        x = r2 * np.cos(t)
        y = r2 * np.sin(t)
        return x, y

    def kiSolver(self):
        R = self.R

        t4 = np.linspace(self.i4,self.i3,self.R)
        L4 = np.trapz(self.rh(t4,self.h2),t4)

        r2 = lambda i: self.r2(i)
        i4 = self.i4
        oi = self.oi

        def fkiSolver(i):
            for val in i:
                kt = np.linspace(i4, val, R)
                L2 = np.trapz(r2(kt), kt*np.sin(oi))
                return abs(L4 - L2)

        ki = minimize(fkiSolver, self.i3).x[0]
        kii = ki * np.sin(oi)
        return ki, kii

    def Pline(self,i,p1,p2):
        a = (p2[1] - p1[1]) / (p2[0] - p1[0])
        b = p1[0]
        c = p1[1]
        return (c - a * b) / (np.sin(i) - a * np.cos(i))

    def Plots(self):
        x1, y1 = self.r1_line()
        x2, y2 = self.r2_line()

        x5 = np.linspace(x1[0],x2[0])
        y5 = np.linspace(y1[0],y2[0])

        x6 = np.linspace(x1[-1],x2[-1])
        y6 = np.linspace(y1[-1],y2[-1])

        p1 = [x1[0],y1[0]]
        p2 = [x2[0],y2[0]]
        t7 = np.linspace(self.ii2, self.ii4, self.R)
        x7 = self.Pline(t7,p1,p2)*np.cos(t7)
        y7 = self.Pline(t7,p1,p2)*np.sin(t7)

        t3 = np.linspace(self.i2,self.i1,self.R)
        x3 = self.r3 * np.cos(t3)
        y3 = self.r3 * np.sin(t3)

        t4 = np.linspace(self.i4, self.i3, self.R)
        x4 = self.rh(t4, self.h2) * np.cos(t4)
        y4 = self.rh(t4, self.h2) * np.sin(t4)

        plots = [[x1,x2,x3,x4,x5,x6,x7],[y1,y2,y3,y4,y5,y6,y7]]

        i = 0
        while i < len(plots[0]):
            plt.plot(plots[0][i], plots[1][i])
            i += 1

        plt.title("")
        plt.grid(True)
        lim = 2
        plt.xlim(4, 10)
        plt.ylim(-lim, lim)
        plt.show()

    def Plot3D(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        dz = 5
        zi = 0
        while zi < len(self.z):
            t = np.linspace(self.ii_if[zi][0],self.ii_if[zi][1],self.R)
            x = self.rh(t,self.z[zi])*np.cos(t)-self.a
            y = self.rh(t, self.z[zi]) * np.sin(t)
            z = np.array([self.z[zi]]*self.R)
            x_back = np.array([-x[0],x[0]])
            y_back = np.array([self.b2]*2)
            z_back = np.array([self.z[zi]]*2)
            ax.plot(x, y, z)
            ax.plot(-x, y, z)
            ax.plot(x_back, y_back, z_back)
            zi += dz

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        lim = 1
        xadd = 0
        yadd = 1+self.b2
        zadd = 1
        ax.set_xlim([-lim+xadd,lim+xadd])
        ax.set_ylim([-lim+yadd, lim+yadd])
        ax.set_zlim([-lim+zadd, lim+zadd])
        plt.show()

    def ACalculator(self):
        t1 = np.linspace(self.ii2, self.ii1, self.R)

        A = np.trapz(0.5 * (self.r2(t1) * self.r2(t1) - self.r1 * self.r1), t1)

        tA = np.linspace(self.ii1, self.kii, self.R)

        p1 = [self.r1*np.cos(self.ii1),self.r1*np.sin(self.ii1)]
        p2 = [self.r2(self.kii)*np.cos(self.kii), self.r2(self.kii)*np.sin(self.kii)]

        A += np.trapz(0.5 * (self.r2(tA) * self.r2(tA) - self.Pline(tA,p1,p2) * self.Pline(tA,p1,p2)), tA)
        self.A = A

    def Get3D(self):
        zi = 0
        dr = 0.004
        d_m = 700
        dz = self.h2 / self.R
        i_i = self.i2
        self.z = []
        self.ii_if = []
        self.LA = []

        self.Cy = 0
        self.Cz = 0
        self.A = 0
        self.P = 0
        Ltot = 0

        while zi < self.h2:
            i_i = minimize(lambda i: abs(self.rh(i, zi) * np.sin(i) - self.b2), i_i).x[0]
            i_f = np.arctan((self.b1 + zi * np.tan(self.of)) / self.a)

            t = np.linspace(i_i, i_f, self.R)

            L = np.trapz(self.rh(t, zi), t)
            Ltot += L
            self.A += np.trapz(self.rh(t,zi)*dz/np.cos(self.op(t)),t)

            x = self.rh(t,zi)*np.cos(t) - self.a
            y = self.rh(t,zi)*np.sin(t)
            A = np.trapz(x,y)

            self.Cy += np.trapz(self.rh(t, zi) * self.rh(t, zi) * np.sin(t), t)
            self.Cz += L*zi

            self.z += [zi]
            self.ii_if += [[i_i, i_f]]
            self.LA += [[L,A]]

            zi += dz

        self.Cy /= Ltot
        self.Cz /= Ltot
        self.P = self.A*dr*d_m

    def Impulsion(self, kg):
        zi = 0
        dz = self.h2/self.R
        da = 997
        Impulsion = 0
        while Impulsion < kg:
            Impulsion += self.LA[zi][1]*dz*da
            zi += 1
        print(f"Linha d'água aos {kg}Kg - {self.z[zi]}m")

    def Save(self):
        save_dict = {
            "r3": self.r3,
            "a": self.a,
            "h2": self.h2,
            "R": self.R,
            "b0": self.b0,
            "b1": self.b1,
            "b2": self.b2,
            "oi": self.oi,
            "of": self.of,
            "ofp": self.ofp,
            "i1": self.i1,
            "i2": self.i2,
            "i3": self.i3,
            "i4": self.i4,
            "l4": self.l4,
            "r1": self.r1,
            "ii1": self.ii1,
            "ii2": self.ii2,
            "ii4": self.ii4,
            "ki": self.ki,
            "kii": self.kii,
            "z": self.z,
            "ii_if": self.ii_if,
            "LA": self.LA,
            "Cy": self.Cy,
            "Cz": self.Cz,
            "A": self.A,
            "P": self.P
        }
        with open('Laterais.json', 'w') as f:
            json.dump(save_dict, f)
        print("Data Saved")

    def Load(self):
        with open('Laterais.json', 'r') as f:
            loaded_data = json.load(f)
        data = loaded_data
        self.r3 = data["r3"]
        self.a = data["a"]
        self.h2 = data["h2"]
        self.R = data["R"]
        self.b0 = data["b0"]
        self.b1 = data["b1"]
        self.b2 = data["b2"]
        self.oi = data["oi"]
        self.of = data["of"]
        self.ofp = data["ofp"]
        self.i1 = data["i1"]
        self.i2 = data["i2"]
        self.i3 = data["i3"]
        self.i4 = data["i4"]
        self.l4 = data["l4"]
        self.r1 = data["r1"]
        self.ii1 = data["ii1"]
        self.ii2 = data["ii2"]
        self.ii4 = data["ii4"]
        self.ki = data["ki"]
        self.kii = data["kii"]
        self.z = data["z"]
        self.ii_if = data["ii_if"]
        self.LA = data["LA"]
        self.Cy = data["Cy"]
        self.Cz = data["Cz"]
        self.A = data["A"]
        self.P = data["P"]
        print("Data Loaded")


l1 = 0.6
l2 = 0.5
h1 = 2

h2 = 0.4
oi = 30 * np.pi/180
of = 45 * np.pi/180

BS = BaseSolver(l1,l2,h1)
LS = LateralSolver()

if 0:
    LS.Build(BS, h2, oi, of)
    LS.Get3D()
    LS.Save()
else:
    LS.Load()
#LS.Impulsion(0.05,yp)
#LS.Plots()