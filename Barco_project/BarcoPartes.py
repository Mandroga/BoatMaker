import numpy as np
from Base import BaseSolver
from Laterais2 import LateralSolver
from Base_de_Mastro import BaseMastroSolver
from Pecas_verticais import VerticaisSolver
import matplotlib.pyplot as plt

l1 = 0.5
l2 = 0.4
h1 = 2

h2 = 0.4
oi = 30 * (np.pi / 180)
of = 45 * (np.pi / 180)




BS = BaseSolver(l1, l2, h1)

LS = LateralSolver(BS, h2, oi, of)

ym = (LS.b0-LS.b2)*(19/26)+LS.b2 # Racio do Laser
rm = 0.03
u0 = 0.05

yp = [0, BS.b1/2]
hyf = 0.05

BMS = BaseMastroSolver(rm,ym, u0, LS)
VS = VerticaisSolver(oi, l2, h2, ym, rm, u0, LS, BMS)

#LS.Impulsion(hyf, yp)
Atotal = BS.A + 2*LS.A + BMS.A + VS.A
print(f"Comprimento maximo = {LS.b0-LS.b2}, Largura maxima = {2*max(LS.x4)}")
print(f"A Total = {Atotal}")


if 0:
    ds = 0.6
    i = 0
    while i<len(BS.plot[0]):
        plt.plot(BS.plot[0][i], BS.plot[1][i])
        i += 1

    i = 0
    while i<len(LS.plot[0])-2:
        plt.plot(LS.plot[0][i]+ds, LS.plot[1][i])
        plt.plot(LS.plot[0][i]+2*ds, LS.plot[1][i])
        i += 1



    i = 0
    while i<len(BMS.plot[0]):
        plt.plot(BMS.plot[0][i]+ds, BMS.plot[1][i])
        i += 1

    i = 0
    while i<len(VS.plot[0]):
        plt.plot(VS.plot[0][i]+ds, VS.plot[1][i])
        i += 1


plt.title("")
ax = plt.gca()
ax.set_aspect('equal')
plt.grid(True)
lim = 4
plt.xlim(-2, 4)
plt.ylim(-2, 4)
plt.show()