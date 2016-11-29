
# coding: utf-8

# In[1]:

# импортируем всё
import numpy as np
import math
from tempreture import tfromT
from nTOp import lambda_n__p, lambda_p__n
from constants import *
import scipy
from scipy import integrate
from scipy.misc import derivative

from rates import p_n__g_d
from tempreture import tfromT

import matplotlib as mpl
import matplotlib.pyplot as plt


# Начальные значения:
# Половина протонов, половина нейтронов, остального нет совсем

# In[2]:

X_0 = np.array([0.5, 0.5, 0.0]) #n, n, d


# Пишем систему диффуров. Полагаем сейчас, что дейтерий лишь образуется. Его распад всё равно не может существенно повлиять на количество протонов и нейтронов
# $$ \left\{
# \begin{aligned}
#     \frac{d X_n}{d T} &= -X_n\lambda_{n\to p}(T)+X_p \lambda_{p \to n} - X_n X_p [np]\\
#     \frac{d X_p}{dT} &= X_n\lambda_{n\to p}(T)-X_p \lambda_{p \to n} - X_n X_p [np]\\
#     \frac{d X_d}{dT} &=2X_n X_p [np]\\
# \end{aligned}
# \right. $$
# 
# 
# Замена **T $\to$ -T** из-за передачи отрицательной последовательности (см. ниже)

# In[3]:

# температуру считаем передаём в эВ
def odu(X, T):
    T = -T
    Xn, Xp, Xd = X
    dxn = -Xn*lambda_n__p(T)+Xp*lambda_p__n(T)-Xp*Xn*p_n__g_d(T)
    dxp = +Xn*lambda_n__p(T)-Xp*lambda_p__n(T)-Xp*Xn*p_n__g_d(T)
    dxd = 2*Xp*Xn*p_n__g_d(T)
    return np.array([dxn, dxp, dxd])

def odu_int(T,X):
    # для odeint
    return odu(X,T)

def jac(X, T):
    T = -T
    Xn, Xp, Xd = X
    #dxn=
    J11 = -lambda_n__p(T) - Xp*p_n__g_d(T) #n
    J12 = lambda_p__n(T) - Xn*p_n__g_d(T) #p
    J13 = 0 #d
    
    #dxp=
    J21 = lambda_n__p(T) - Xp*p_n__g_d(T) #n
    J22 = -lambda_p__n(T) - Xn*p_n__g_d(T) #p
    J23 = 0 #d
    
    #dxd=
    J31 = 2*Xp*p_n__g_d(T) #n
    J32 = 2*Xn*p_n__g_d(T) #p
    J33 = 0 #d

    return [
        [J11, J12, J13],
        [J21, J22, J23],
        [J31, J32, J33]
    ]

def jac_i(T, X):
    return jac(X,T)


# Теперь разберёмся с сеткой. У нас есть температура, в ней мы найдём ответ (т.е. dX/dT) и уже конечный преобразуем во время

# In[4]:

Ts = np.logspace(math.log10(10**11), math.log10(9*10**10), num=100)*k_b


# Поскольку нам нужно аналог возрастающей последовательности, а *Ts* - убывает, домножим на **-1** последовательность и учтём это в функциях

# In[5]:

Ts = -Ts


# In[6]:

odes = integrate.ode(odu_int, jac=jac_i)
odes.set_integrator('vode', method="bdf", order=4, nsteps=10000)
odes.set_initial_value(X_0, Ts[0])


# In[7]:

sol_n, sol_p, sol_d, tres = [],[], [], []
tres.append(Ts[0])
sol_n.append(X_0[0])
sol_p.append(X_0[1])
sol_d.append(X_0[2])
i=0
while odes.successful() and odes.t < Ts[-1]:
    dt = Ts[i+1]-Ts[i]
    solu = odes.integrate(odes.t+dt)
    i+=1
    sol_n.append(solu[0])
    sol_p.append(solu[1])
    sol_d.append(solu[2])
    tres.append(odes.t+dt)


# In[8]:

plt.plot(tres, sol_n, 'b', label='n(t)')
plt.plot(tres, sol_p, 'g', label='p(t)')
plt.plot(tres, sol_d, 'r', label='d(t)')


# In[9]:

plt.show()


# In[10]:

get_ipython().system('notify-send "hello"')


# In[11]:

import numpy as np
import math
from tempreture import tfromT
from nTOp import lambda_n__p, lambda_p__n
from constants import *
import scipy
from scipy import integrate
from scipy.misc import derivative
from rates import p_n__g_d

import matplotlib as mpl
import matplotlib.pyplot as plt

_X_n_0 = 0.5
_X_p_0 = 0.5
_X_0 = np.array([0.5, 0.5, 0.0])
# def g(y,x):
    # y0

Ts = np.logspace(math.log10(10**11), math.log10(10**10), num=40)*k_b
Ts = -Ts
def ode_int(X, T):
    T = -T
    X_n, X_p, X_d = X
    lam_np = lambda_n__p(T)
    lam_pn = lambda_p__n(T)
    pnd = p_n__g_d(T)
    if pnd > min(lam_np, lam_pn):
        pnd = min(lam_np, lam_pn)/10**4
    dX_n = (-lam_np*X_n+lam_pn*(X_p)) - X_p*X_n*pnd
    dX_p = (lam_np*X_n-lam_pn*(X_p)) - X_p*X_n*pnd
    dX_d = 2*X_p*X_n*p_n__g_d(T)
    global num
    if num % 1000 == 0:
        print(lam_np, lam_pn, pnd)
    return [dX_n, dX_p, dX_d] #dX_n #
num = 0
def ode_(T, X):
    global num
    if num % 1000 == 0:
        print(num)
    num += 1
    return ode_int(X, T)

def jacob(T,X):
    T = -T
    X_n, X_p, X_d = X
    return [
        [-lambda_n__p(T) - X_p*p_n__g_d(T), +lambda_p__n(T) - X_n*p_n__g_d(T), 0],
        [+lambda_n__p(T) - X_p*p_n__g_d(T), -lambda_p__n(T) - X_n*p_n__g_d(T), 0],
        [2*X_p*p_n__g_d(T), 2*X_p*p_n__g_d(T), 0]
    ]


# In[12]:

odes = integrate.ode(ode_, jac=jacob)
odes.set_integrator('vode', method="bdf", nsteps=3000)
odes.set_initial_value(_X_0, Ts[0])
X_ans = _X_0.reshape((1,-1))
Tres=[Ts[0]]
i = 0
while odes.successful() and odes.t < Ts[-1]:
    dt = Ts[i+1]-Ts[i]
    solu = np.array(list(odes.integrate(odes.t+dt))).reshape((1,-1))
    print(solu)
    i+=1
    X_ans = np.append(X_ans, solu, axis=0)
    Tres.append(odes.t+dt)


# In[23]:

X_ans


# In[14]:

ts = list(map(tfromT, -np.array(Tres)))
# %matplotlib inline
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
# plt.yscale('log')
plt.xlabel(r'\textbf{t} (s)')
plt.ylabel(r'\textbf{X_n}')

plt.plot(ts, X_ans[:,0], 
        linewidth=2.0, label=r'X_n')

plt.plot(ts, X_ans[:,1], 
        linewidth=2.0, label=r'X_p')
plt.plot(ts, X_ans[:,1], 
        linewidth=2.0, label=r'X_d')

plt.legend()
#plt.gca().invert_xaxis()
plt.show()


# In[15]:

Ts/k_b


# In[ ]:



