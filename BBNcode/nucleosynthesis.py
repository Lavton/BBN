"""
"Влияние тяжёлых изотопов гелия на процессы первичного нуклеосинтеза"
"""
import logging
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d, FUNC:%(funcName)s]# %(levelname)-8s  %(message)s', 
    level=logging.DEBUG,
    filename=u'mylog.log')

import numpy as np
import time
import math
import pickle
from tempreture import tfromT, Tfromt, derriviate_T_from_t
from nTOp import lambda_n__p, lambda_p__n
import scipy
from scipy import integrate
from scipy.misc import derivative
import constants
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import elements.register as elements
import datetime
import os
now_title = datetime.datetime.now().isoformat()

start_time = time.time()

# начальные массовые доли элементов берём из elements
X_0 = np.array(elements.X_0)

# обезразмеренный диапазон температур
grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
grid2 = grid
Ts = constants.less_tempreture(grid2, units="K")
# переводим в отрицательную шкалу, чтобы Ts[i] > Ts[i-1]
ts = np.array([tfromT(T) for T in Ts])

with open("mylog.log", "wt") as f:
    f.write("")
logging.info(("start", constants.ode_params, ts))

def ode_int(X, t):
    """
    вычисляем следующий шаг для интегрирования
    X - значения массовых долей элементов на предыдущем шаге, T - температура

    тут выполняем технические части. Логика - в @see elements.registrator.sode_int
    вовзвращаем изменения элементов
    """
    dx = elements.registrator.sode_int(X=X, T=Tfromt(t))
    dX = np.array(dx)
    return dX

num = 0
def ode_(t, X):
    """
    то же, что @see ode_int, только входные параметры T и X поменены местами
    """
    global num
    if num % 100 == 0:
        logging.debug(num)
    num += 1
    return ode_int(X, t)


def jacob(t, X):
    """
    вычисляем Якобиан уравнения. Логика спрятана в @elements.registrator.jacob
    """
    j = np.array(elements.registrator.jacob(X, Tfromt(t))) # * derriviate_T_from_t(T)

    return j


class TechnicalCalcExitException(Exception):
    def __init__(self, i, Tres, X_ans):
        super(Exception, self).__init__()
        self.i = i
        self.Tres = Tres
        self.X_ans = X_ans


def technical_stop_cache(i, Tres, X_ans):

    with open("smart_cache.pickle", "wb") as f:
        pickle.dump((i, X_ans, Tres, constants.ode_params), f)
    if os.path.isfile("exit_now"):
        os.remove("exit_now")
        raise TechnicalCalcExitException(i, Tres, X_ans)


# инициируем программу для решение дифура
def iter_process(X_0, T0, Ts, i, X_ans, Tres):
    # выполняем шаги
    for param_set in constants.ode_params:
        if Tres[-1] >= param_set[0]:
            if "rtoi" in param_set[1]:
                constants.def_params["rtoi"] = param_set[1]["rtoi"]
            if "max_step" in param_set[1]:
                constants.def_params["max_step"] = param_set[1]["max_step"]
            if "min_step" in param_set[1]:
                constants.def_params["min_step"] = param_set[1]["min_step"]
            if "atoi" in param_set[1]:
                constants.def_params["atoi"] = param_set[1]["atoi"]

    last_step = Tres[-1]
    logging.debug("new start iter_process with i = {}, t = {}, X = {}".format(i, T0, X_0))
    odes = integrate.ode(ode_, jac=jacob)
    odes.set_integrator('vode', method="bdf", with_jacobian=True, nsteps=8000, 
        min_step=constants.def_params["min_step"],
        rtol=constants.def_params["rtoi"],
        max_step=constants.def_params["max_step"],
        atol=constants.def_params["atoi"]
        )
    odes.set_initial_value(X_0, T0)
    while odes.successful() and odes.t < Ts[-2]:
        technical_stop_cache(i, Tres, X_ans)
        if abs(sum(X_ans[-1])-1.0) >= 1e-6:
            logging.warning(("X error to big!", abs(sum(X_ans[-1])-1.0)))
        for param_set in constants.ode_params:
            if Tres[-1] >= param_set[0]:
                if param_set[0] >= last_step:
                    logging.debug(("exit on param_set", param_set))
                    return (i, X_ans, Tres)
        for element in elements.registrator.elements:
            # во избежание численных ошибок, концентрация элементов изначально считается из
            # закона равновесия, и лишь потом входит в полноценный диффур
            if element.equilibrium:
                if Tres[-1] >= element.tr_t:
                    if element.tr_t >= last_step:
                        logging.info(("start ode of element", element.str_view))
                        return (i, X_ans, Tres)

        dt = Ts[i+1]-Ts[i]
        step_solu = odes.integrate(odes.t+dt)
        solu = np.array(list(step_solu)).reshape((1,-1))
        Tres.append(odes.t+dt)
        for element in elements.registrator.elements:
            # во избежание численных ошибок, концентрация элементов изначально считается из
            # закона равновесия, и лишь потом входит в полноценный диффур
            if element.equilibrium:
                if Tres[-2] >= last_step >= element.tr_t:
                    logging.debug("noteq")
                    pass
                else:
                    logging.debug(("doeq", element.str_view))
                    solu = element.equilibrium(solu, Tfromt(odes.t))
        i+=1
        logging.info("i = {}/{}, t = {}, last X = {}".format(i, len(Ts), Tres[-1], X_ans[-1]))
        X_ans = np.append(X_ans, solu, axis=0)
    return(-1, X_ans, Tres)

X_ans = X_0.reshape((1,-1))
i = 0
Tres = [ts[i]]

def start_from_cache(i, X_ans, Tres):
    want_t = 1e30 
    if len(sys.argv) > 1:
        want_t = eval(sys.argv[1])  # с какого времени начинать
    if constants.smart_caching and os.path.isfile("smart_cache.pickle"):
        with open("smart_cache.pickle", "rb") as f:
            t_i, t_X_ans, t_Tres, t_ode_params = pickle.load(f)
        for j in range(min(len(t_ode_params), len(constants.ode_params))):
            if t_ode_params[i] != constants.ode_params[i]:
                break
        for ii in range(t_i):
            if t_Tres[ii] >= constants.ode_params[j][0]:
                break
            if t_Tres[ii]>= want_t:
                break
        logging.info("use smart cache. Start from {}".format(ii))
        return (ii, t_X_ans[:ii+1], t_Tres[:ii+1])
    else:
        return (i, X_ans, Tres)
        
i, X_ans, Tres = start_from_cache(i, X_ans, Tres)

try:
    while i != -1:
        X_0 = X_ans[-1]
        i, X_ans, Tres = iter_process(X_0, ts[i], ts, i, X_ans, Tres)
except TechnicalCalcExitException as e:
    i, Tres, X_ans = e.i, e.Tres, e.X_ans
    logging.warning("interrupt by user")


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
logging.info(("TIME WORKS", (time.time() - start_time)/60))
# plt.title(now_title)
with open("Output.pickle", "wb") as f:
    pickle.dump((X_ans, Tres), f)

sys.stdout.flush()
plt.show()
