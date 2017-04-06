"""
"Влияние тяжёлых изотопов гелия на процессы первичного нуклеосинтеза"
"""

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
    dx = elements.registrator.sode_int(X=X, T=T)
    dX = np.array(dx) * derriviate_T_from_t(T)
    T = -T
    return dX

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


odes = integrate.ode(ode_, jac=jacob)
odes.set_integrator('vode', method="bdf", nsteps=2000)
odes.set_initial_value(X_0, Ts[0])
X_ans = X_0.reshape((1,-1))
Tres=[Ts[0]]
i = 0
while odes.successful() and odes.t < Ts[-1]:
    dt = Ts[i+1]-Ts[i]
    solu = np.array(list(odes.integrate(odes.t+dt))).reshape((1,-1))
    odes._y[-1] = 1e-12
    i+=1
    Tres.append(odes.t+dt)
    for element in elements.registrator.elements:
        if element.equilibrium:
            if -Tres[-1] > element.tr_T:
                solu = element.equilibrium(solu, Tres[-1])
            else:
                if not element.is_ode_state:
                    element.is_ode_state = True
                    solu = element.equilibrium(solu, Tres[-1])
                    odes = integrate.ode(ode_, jac=jacob)
                    odes.set_integrator('vode', method="bdf", nsteps=800)
                    odes.set_initial_value(solu[0], Tres[-1])
                    print("Here", solu[0])


    print(solu, "i = {}/{}".format(i, len(Ts)))
    X_ans = np.append(X_ans, solu, axis=0)

ts = [constants.to_norm_time(t) for t in  map(tfromT, -np.array(Tres))]
# %matplotlib inline
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\textbf{t} (s)$')
plt.xlim([1e-2,1e3])
ylabel = r"\textbf{X}"
plt.ylabel(ylabel)
plt.ylim([1e-14, 3])
elements.registrator.calc_plot(plt, ts, X_ans)
from elements.H_2 import H_2
plt.axvline(x=constants.to_norm_time(tfromT(H_2.tr_T)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{tempreture} (MeV)')
plt.ylabel(r'\textbf{\lambda}')

# plt.legend()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
#plt.gca().invert_xaxis()
plt.show()