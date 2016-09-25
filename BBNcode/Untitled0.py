
# coding: utf-8

# In[19]:

# %load __init__.py
import numpy as np
import math
from tempreture import tfromT
from nTOp import lambda_n__p, lambda_p__n
from constants import *
from scipy import integrate
from scipy.misc import derivative
from rates import p_n__g_d

import matplotlib as mpl
import matplotlib.pyplot as plt

_X_n_0 = 0.5
_X_p_0 = 0.5
_X_d_0 = 0.0
X_0 = np.array([_X_n_0, _X_p_0, _X_d_0])

# def g(y,x):
    # y0

Ts = np.logspace(math.log10(10**11), math.log10(10**8), num=20)*k_b

func = lambda X_n, T: (-lambda_n__p(T)*X_n+lambda_p__n(T)*(1-X_n))*(
    1/derivative(tfromT, T))
def func(X, T):
	return np.array([(-X[0]*lambda_n__p(T)+lambda_p__n(T)*X[1])*(
    1/derivative(tfromT, T)),
    (+X[0]*lambda_n__p(T)-lambda_p__n(T)*X[1])*(
    1/derivative(tfromT, T)),
    (X[0]*X[1]*p_n__g_d(T))*(
    1/derivative(tfromT, T))])



# In[20]:

X = integrate.odeint(func, X_0, Ts)



# In[30]:

X_n


# In[22]:

X_p = 1 - X_n


# In[31]:

ts = list(map(tfromT, Ts))
ts


# In[34]:

get_ipython().magic('matplotlib inline')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{t} (s)')
plt.ylabel(r'\textbf{X_n}')

plt.plot(ts, X_n, 
        linewidth=2.0, label=r'X_n')

plt.plot(ts, X_p, 
        linewidth=2.0, label=r'X_p')

plt.legend()
#plt.gca().invert_xaxis()
plt.show()


# In[ ]:



