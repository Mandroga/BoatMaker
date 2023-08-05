import numpy as np
from Laterais2 import LateralSolver
from Base import BaseSolver
import matplotlib.pyplot as plt
from Tools import ComputationalMath as CM

l1 = 0.4
l2 = 0.3
h1 = 2

h2 = 0.3
oi = 30 * np.pi/180
of = 45 * np.pi/180

BS = BaseSolver(l1,l2,h1)

LS = LateralSolver()
LateralSolver.Build(LS,BS, h2, oi, of)

yp = [0, BS.b1/2]
hyf = 0.05
LS.Impulsion(hyf, yp)

plt.plot(LS.x4,LS.y4)
plt.plot(-LS.x4,LS.y4)
for val in yp:
    plt.plot([-LS.l4/2,LS.l4/2],[val]*2)
plt.plot([-LS.l4/2,LS.l4/2],[LS.b2]*2)

plt.title("")
ax = plt.gca()
#ax.set_aspect('equal')
plt.grid(True)
lim = 2
plt.xlim(-lim,lim)
plt.ylim(-lim, lim)
plt.show()