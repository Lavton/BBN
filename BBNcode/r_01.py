# импортируем всё
import numpy as np
import math
from tempreture import tfromT, derriviate_T_from_t
from nTOp import lambda_n__p, lambda_p__n
from constants import *
import constants
import scipy
from scipy import integrate
from scipy.misc import derivative

from rates import p_n__g_d
from tempreture import tfromT

import matplotlib as mpl
import matplotlib.pyplot as plt

from collections import namedtuple

Elements = namedtuple('Elements', ['n', 'p'])

X_0 = np.array([0.5, 0.5])


def odu(X, T):
    # T = -T
    Xn, Xp = X
    dxn = (-Xn*lambda_n__p(T)+Xp*lambda_p__n(T)) * derriviate_T_from_t(T)
    dxp = +Xn*lambda_n__p(T)-Xp*lambda_p__n(T) * derriviate_T_from_t(T)
    return np.array([dxn, dxp])

def odu_int(T,X):
    # для odeint
    return odu(X,T)

def jac(X, T):
    # T = -T
    Xn, Xp, Xd = X
    #dxn=
    J11 = (-lambda_n__p(T))  * derriviate_T_from_t(T) #n
    J12 = (lambda_p__n(T)) * derriviate_T_from_t(T) #p
    # J13 = 0 #d
    
    #dxp=
    J21 = (lambda_n__p(T)) * derriviate_T_from_t(T) #n
    J22 = (-lambda_p__n(T)) * derriviate_T_from_t(T) #p
    # J23 = 0 #d
    
    #dxd=
    # J31 = 2*Xp*p_n__g_d(T) #n
    # J32 = 2*Xn*p_n__g_d(T) #p
    # J33 = 0 #d

    return [
        [J11, J12],
        [J21, J22],
        # [J31, J32, J33]
    ]

def jac_i(T, X):
    return jac(X,T)

Ts = constants.less_tempreture(np.logspace(math.log10(10**11), math.log10(10**8), num=100), units="K")
odes = integrate.ode(odu_int, jac=jac_i)
odes.set_integrator('vode', method="bdf", order=4, nsteps=300)
odes.set_initial_value(X_0, Ts[0])

sol_n, sol_p, tres = [],[], []
tres.append(Ts[0])
sol_n.append(X_0[0])
sol_p.append(X_0[1])
# sol_d.append(X_0[2])
i=0
while odes.successful() and odes.t < Ts[-1]:
    dt = Ts[i+1]-Ts[i]
    solu = odes.integrate(odes.t+dt)
    i+=1
    sol_n.append(solu[0])
    sol_p.append(solu[1])
    # sol_d.append(solu[2])
    tres.append(odes.t+dt)

plt.plot(tres, sol_n, 'b', label='n(t)')
plt.plot(tres, sol_p, 'g', label='p(t)')
# plt.plot(tres, sol_d, 'r', label='d(t)')
plt.show()