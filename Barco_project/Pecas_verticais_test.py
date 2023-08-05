import numpy as np
from Laterais2 import LateralSolver
from Base import BaseSolver
from Base_de_Mastro import BaseMastroSolver
from Pecas_verticais import VerticaisSolver
import matplotlib.pyplot as plt

l1 = 0.4
l2 = 0.3
h1 = 2

BS = BaseSolver(l1, l2, h1)
b2 = BS.b2

h2 = 0.3
oi = 30 * (np.pi / 180)
of = 45 * (np.pi / 180)

LS = LateralSolver(BS,h2,oi,of,10**-3,100)
b0 = LS.b0

ym = (b0-b2)*(19/26)+b2
rm = 0.03
u0 = 0.05
BMS = BaseMastroSolver(rm,ym,u0,LS)
VS = VerticaisSolver(oi, l2, h2, ym, rm, u0, LS, BMS)

i = 0
while i<len(VS.plot[0]):
    plt.plot(VS.plot[0][i], VS.plot[1][i])
    i += 1

plt.title("")
plt.grid(True)
lim = 5
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.show()
