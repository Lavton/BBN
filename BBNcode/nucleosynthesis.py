import numpy as np
import math
from tempreture import tfromT, derriviate_T_from_t
from nTOp import lambda_n__p, lambda_p__n
import scipy
from scipy import integrate
from scipy.misc import derivative
import constants

import matplotlib as mpl
import matplotlib.pyplot as plt
import elements.register as elements


X_0 = np.array(elements.X_0)

print(X_0)

Ts = constants.less_tempreture(np.logspace(math.log10(10**11), math.log10(10**9), num=40), units="K")
Ts = -Ts

def ode_int(X, T):
    T = -T
    # X_n, X_p, X_d = X
    # dX_n = (-lambda_n__p(T)*X_n+lambda_p__n(T)*(X_p)) * derriviate_T_from_t(T)
    # dX_p = (lambda_n__p(T)*X_n-lambda_p__n(T)*(X_p))  * derriviate_T_from_t(T)
    # dX_d = 0
    dX = np.array(elements.registrator.ode_int(X, T)) * derriviate_T_from_t(T)
    T = -T
    return dX
    # return [dX_n, dX_p, dX_d]
    # 

num = 0
def ode_(T, X):
    global num
    if num % 100 == 0:
        print(num)
    num += 1
    return ode_int(X, T)


def jacob(T,X):
    T = -T
    j = np.array(elements.registrator.jacob(X, T)) * derriviate_T_from_t(T)

    T = -T
    return j