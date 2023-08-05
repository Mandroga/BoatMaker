import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from Tools import ComputationalMath as CM
from Base import BaseSolver
import time
R = 100
x = np.array([1]*R)
y = np.linspace(0,10,R)

A = np.trapz(x,y)
y_bar = np.trapz(x*y,y)/A

print(f"A = {A}, y_bar = {y_bar}")
plt.plot(x,y)
plt.grid(True)
plt.show()