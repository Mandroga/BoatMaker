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
LS = LateralSolver(BS, h2, oi, of)

#plt.plot(LS.t2,LS.k)
#a, b, c= LS.fit_params
#plt.plot(LS.kt,LS.kfit(LS.kt, a, b, c))

i = 0
while i<len(BS.plot[0]):
    #plt.plot(BS.plot[0][i], BS.plot[1][i])
    i += 1

i = 0
while i<len(LS.plot[0]):
    plt.plot(LS.plot[0][i], LS.plot[1][i])
    i += 1

i = 0
while i<len(LS.plot[0]):
    #plt.plot(LS.plot[0][i]+0.5, LS.plot[1][i])
    i += 1

plt.title("")
ax = plt.gca()
#ax.set_aspect('equal')
plt.grid(True)
lim = 2
plt.xlim(-lim,lim)
plt.ylim(-lim, lim)
plt.show()