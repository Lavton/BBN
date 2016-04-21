import numpy as np
import math
from tempreture import tfromT
from nTOp import lambda_n__p, lambda_p__n
from constants import *
from scipy import integrate
from scipy.misc import derivative

import matplotlib as mpl
import matplotlib.pyplot as plt

_X_n_0 = 0.5
_X_p_0 = 0.5

# def g(y,x):
    # y0

Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=20)

func = lambda X_n, T: (-lambda_n__p(T, units="K")*X_n+lambda_p__n(T, units="K")*(1-X_n))*(derivative(lambda t: tfromT(t, units="K"), T))
X_n = integrate.odeint(func, _X_n_0, Ts)
X_p = 1 - X_n

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{tempreture} (K)')
plt.ylabel(r'\textbf{X_n}')

plt.plot(Ts, X_n, 
        linewidth=2.0, label=r'X_n')

plt.plot(Ts, X_p, 
        linewidth=2.0, label=r'X_p')

plt.legend()
plt.gca().invert_xaxis()
plt.show()