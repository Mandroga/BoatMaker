import numpy as np
from Base import BaseSolver
import matplotlib.pyplot as plt

l1 = 0.4
l2 = 0.3
h1 = 2

BS = BaseSolver(l1,l2,h1)

i = 0
while i<len(BS.plot[0]):
    plt.plot(BS.plot[0][i], BS.plot[1][i])
    i += 1

plt.title("")
plt.grid(True)
lim = BS.b1+0.1
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.show()
