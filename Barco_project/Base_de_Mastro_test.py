import numpy as np
from Base import BaseSolver
from Laterais2 import LateralSolver
from Base_de_Mastro import BaseMastroSolver
import matplotlib.pyplot as plt

l1 = 0.4
l2 = 0.3
h1 = 2

BS = BaseSolver(l1, l2, h1)
b2 = BS.b2

h2 = 0.3
oi = 30 * (np.pi / 180)
of = 45 * (np.pi / 180)

LS = LateralSolver(BS,h2,oi,of,10**-3,1000)
b0 = LS.b0

rm = 0.03
ym = (b0-b2)*(19/26)+b2
u0 = 0.05

BMS = BaseMastroSolver(rm,ym,u0,LS)

i = 0
while i<len(BMS.plot[0]):
    plt.plot(BMS.plot[0][i], BMS.plot[1][i])
    i += 1

plt.title("")
plt.grid(True)
lim = 5
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.show()