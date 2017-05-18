"""
"Влияние тяжёлых изотопов гелия на процессы первичного нуклеосинтеза"
"""

import numpy as np
import time
import math
from tempreture import tfromT, Tfromt, derriviate_T_from_t
from nTOp import lambda_n__p, lambda_p__n
import scipy
from scipy import integrate
from scipy.misc import derivative
import constants

import matplotlib as mpl
import matplotlib.pyplot as plt
import elements.register as elements
import datetime
now_title = datetime.datetime.now().isoformat()

start_time = time.time()
# начальные массовые доли элементов берём из elements
X_0 = np.array(elements.X_0)

print(X_0)

# обезразмеренный диапазон температур
grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
grid2 = grid
# grid2 = np.array(sorted(list(set(list(np.logspace(math.log10(grid[100]), math.log10(10**7), num=50))+list(grid))), reverse=True))
Ts = constants.less_tempreture(grid2, units="K")
# переводим в отрицательную шкалу, чтобы Ts[i] > Ts[i-1]
ts = np.array([tfromT(T) for T in Ts])
print(ts)
# exit()
# Ts = -Ts

def ode_int(X, t):
    """
    вычисляем следующий шаг для интегрирования
    X - значения массовых долей элементов на предыдущем шаге, T - температура

    тут выполняем технические части. Логика - в @see elements.registrator.sode_int
    вовзвращаем изменения элементов
    """
    # T = -T
    dx = elements.registrator.sode_int(X=X, T=Tfromt(t))
    dX = np.array(dx) # * derriviate_T_from_t(T)
    # T = -T
    return dX

num = 0
def ode_(t, X):
    """
    то же, что @see ode_int, только входные параметры T и X поменены местами
    """
    global num
    if num % 100 == 0:
        print(num)
    num += 1
    return ode_int(X, t)


def jacob(t, X):
    """
    вычисляем Якобиан уравнения. Логика спрятана в @elements.registrator.jacob
    """
    # T = -T
    j = np.array(elements.registrator.jacob(X, Tfromt(t))) # * derriviate_T_from_t(T)

    # T = -T
    return j

ode_params = [
    [-1.0, {
    "rtoi": 1e-6,
    "max_step": 1.0,
    "min_step": 1e-10,
    }], 
    [0.004, {
    "rtoi": 1e-7,
    }],
    [0.014, {
    "rtoi": 1e-11,
    "max_step": 0.0005
    ,
    }],
    [0.02, {
    "atoi": 1e-9,
    "rtoi": 1e-11,
    "max_step": 0.0001,
    "min_step": 0.0
    }],
    [0.038, {
    "atoi": 1e-12,
    "max_step": 0.002,
    }],
    [0.052, {
    "rtoi": 1e-6,
    "max_step": 0.2,
    "min_step": 0.0
    }],
    [0.055,{
    "rtoi": 1e-9,
    "max_step": 0.002,
    "atoi": 1e-8,
    }],
    [0.080,{
    "rtoi": 1e-7,
    # "min_step": 1e-8,
    "max_step": 0.03,
    }],
    [0.2, {
    "rtoi": 1e-8,
    "atoi": 1e-10,
    }],
    [0.31, {
    "rtoi": 1e-11,
    "atoi": 1e-12,
    "max_step": 0.003
    }],
    [10, {
    "rtoi": 1e-9,
    "max_step": 0.0
    }],
    [100, {
    "max_step": 0.0,
    }]
]

# инициируем программу для решение дифура
def iter_process(X_0, T0, Ts, i, X_ans, Tres):
    # выполняем шаги
    odes = integrate.ode(ode_, jac=jacob)
    rtoi = 1e-6
    atoi = 1e-12
    max_step = 1.0
    min_step = 1e-10
    last_step = 0.0
    for param_set in ode_params:
        if Tres[-1] >= param_set[0]:
            if "rtoi" in param_set[1]:
                rtoi = param_set[1]["rtoi"]
            if "max_step" in param_set[1]:
                max_step = param_set[1]["max_step"]
            if "min_step" in param_set[1]:
                min_step = param_set[1]["min_step"]
            if "atoi" in param_set[1]:
                atoi = param_set[1]["atoi"]
    last_step = Tres[-1]
    print("start", Tres[-1])
    odes.set_integrator('vode', method="bdf", with_jacobian=True, nsteps=8000, 
        min_step=min_step, 
        rtol=rtoi, 
        max_step=max_step,
        atol=atoi
        )
    odes.set_initial_value(X_0, T0)
    while odes.successful() and odes.t < Ts[-2]:
        for param_set in ode_params:
            if Tres[-1] >= param_set[0]:
                if param_set[0] >= last_step:
                    return (i, X_ans, Tres)
        for element in elements.registrator.elements:
            # во избежание численных ошибок, концентрация элементов изначально считается из
            # закона равновесия, и лишь потом входит в полноценный диффур
            if element.equilibrium:
                if Tres[-1] >= element.tr_t:
                    if element.tr_t >= last_step:
                        print("element", element.str_view)
                        return (i, X_ans, Tres)

        dt = Ts[i+1]-Ts[i]
        step_solu = odes.integrate(odes.t+dt)
        print(i)
        solu = np.array(list(step_solu)).reshape((1,-1))
        i+=1
        Tres.append(odes.t+dt)
        flag_was_eq = False
        for element in elements.registrator.elements:
            # во избежание численных ошибок, концентрация элементов изначально считается из
            # закона равновесия, и лишь потом входит в полноценный диффур
            if element.equilibrium:
                if Tres[-2] >= element.tr_t >= last_step:
                    pass
                else:
                    solu = element.equilibrium(solu, Tfromt(odes.t))
        print(solu, "i = {}/{}".format(i, len(Ts)), Tres[-1], "raz", sum(solu[0]) - 1)
        X_ans = np.append(X_ans, solu, axis=0)
        if flag_was_eq:
            print("here")
            return (i, X_ans, Tres)
    return(-1, X_ans, Tres)

# odes = integrate.ode(ode_, jac=jacob)
# odes.set_integrator('vode', method="bdf", nsteps=2000)
# odes.set_initial_value(X_0, Ts[0])
X_ans = X_0.reshape((1,-1))
i = 0
# Tres=[Ts[i]]
Tres = [ts[i]]

while i != -1:
    X_0 = X_ans[-1]
    # print(ts[i])
    i, X_ans, Tres = iter_process(X_0, ts[i], ts, i, X_ans, Tres)

# время
# ts = [constants.to_norm_time(t) for t in  map(tfromT, -np.array(Tres))]

# рисуем всё, что можно :)
time.sleep(1)
my_xn = list(X_ans[:,2])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\textbf{t} (s)$')
plt.xlim([1e-4,1e3])
ylabel = r"\textbf{X}"
plt.ylabel(ylabel)
plt.ylim([1e-30, 1000])
elements.registrator.calc_plot(plt, Tres, X_ans)
###################
import elements._xn_modeling_wai
Ts_ = []
Tnus_ = []
tu = []
for (T_, Xn_) in elements._xn_modeling_wai.t_xn.items():
    tu.append((T_, Xn_))

tu = sorted(tu, reverse=True)
tu
for (T_, Xn_) in tu:
    Ts_.append(T_)
    Tnus_.append(Xn_)
    
plt.plot([tfromT(constants.less_tempreture(T, units="K")) for T in Ts_], Tnus_, label="tabular result")


for t in Tres:
  plt.axvline(x=t, linewidth=0.1)  
##################
from elements.H_2 import H_2
plt.axvline(x=H_2.tr_t)
plt.plot(Tres, [H_2.equilibrium([[X_ans[i][0], X_ans[i][1], 0]], Tfromt(Tres[i]))[0][2] for i in range(len(Tres))])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{time} (s)')
plt.ylabel(r'\textbf{\lambda}')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
print("ode_params", ode_params)
print("TIME WORKS", (time.time() - start_time)/60)
plt.title(now_title)
import pickle
with open("Output.pickle", "wb") as f:
    pickle.dump((X_ans, Tres), f)


plt.show()
# time.sleep(1)
# for k in range(len(Tres)):
#     if k:
#         if my_xn[k] <= my_xn[k-1]:
#             print("AAAAAAAAAA")
#     print(Tres[k], my_xn[k])

# print(my_xn)