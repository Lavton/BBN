"""
"Влияние тяжёлых изотопов гелия на процессы первичного нуклеосинтеза"
"""
import logging
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d, FUNC:%(funcName)s]# %(levelname)-8s  %(message)s', 
    level=logging.INFO)
# ,
    # filename=u'mylog.log')

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
import csv
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import register
import datetime
import tempreture
import os
now_title = datetime.datetime.now().isoformat()

start_time = time.time()

# начальные массовые доли элементов берём из register
X_0 = np.array(register.X_0)

# обезразмеренный диапазон температур
grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
grid2 = grid
Ts = constants.less_tempreture(grid2, units="K")
# переводим в отрицательную шкалу, чтобы Ts[i] > Ts[i-1]
ts = np.array([tfromT(T) for T in Ts])

with open("mylog.log", "wt") as f:
    f.write("")

def ode_int(X, t):
    """
    вычисляем следующий шаг для интегрирования
    X - значения массовых долей элементов на предыдущем шаге, T - температура

    тут выполняем технические части. Логика - в @see register.registrator.sode_int
    вовзвращаем изменения элементов
    """
    dx = register.registrator.sode_int(X=X, T=Tfromt(t))
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
    вычисляем Якобиан уравнения. Логика спрятана в @register.registrator.jacob
    """
    j = np.array(register.registrator.jacob(X, Tfromt(t))) # * derriviate_T_from_t(T)

    return j


class TechnicalCalcExitException(Exception):
    def __init__(self, i, Tres, X_ans):
        super(Exception, self).__init__()
        self.i = i
        self.Tres = Tres
        self.X_ans = X_ans


def technical_stop_cache(i, Tres, X_ans):
    with open("smart_cache.pickle", "wb") as f:
        pickle.dump((i, X_ans, Tres, constants.ode_params, [el.names for el in register.registrator.elements]), f)
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
            logging.error(("X error to big!", abs(sum(X_ans[-1])-1.0)))
        if len(X_ans) > 3:
            for kk in range(len(X_ans[-1])):
                if X_ans[-2][kk] == 0:
                    continue
                ma = max(abs(X_ans[-1][kk]), abs(X_ans[-2][kk]))
                mi = min(abs(X_ans[-1][kk]), abs(X_ans[-2][kk]))
                if 20 < ma/mi:
                    logging.error(("dX to big!", register.registrator.elements[kk].str_view, X_ans[-1][kk]/X_ans[-2][kk]))
        if len(list(filter(lambda l: l<0, X_ans[-1]))):
            logging.error(("ONE X LESS ZERRO", X_ans[-1]))
        for param_set in constants.ode_params:
            if Tres[-1] >= param_set[0]:
                if param_set[0] >= last_step:
                    logging.debug(("exit on param_set", param_set))
                    return (i, X_ans, Tres)
        for element in register.registrator.elements:
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
        for element in register.registrator.elements:
            # во избежание численных ошибок, концентрация элементов изначально считается из
            # закона равновесия, и лишь потом входит в полноценный диффур
            if element.equilibrium:
                if Tres[-2] >= last_step >= element.tr_t:
                    logging.debug("noteq")
                    pass
                else:
                    logging.debug("sum X = {}".format(sum(solu[-1])-1))
                    logging.debug(("doeq", element.str_view))
                    solu = element.equilibrium(solu, Tfromt(odes.t))
        i+=1
        logging.info("i = {}/{}, t = {}, T9 = {}, last X = {}".format(
            i, 
            len(Ts), 
            Tres[-1], 
            constants.to_norm_tempreture(tempreture.Tfromt(Tres[-1]), units="T9"),
            X_ans[-1][::-1])
        )
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
            t_i, t_X_ans, t_Tres, t_ode_params, t_elements = pickle.load(f)
        for j in range(min(len(t_ode_params), len(constants.ode_params))):
            if t_ode_params[i] != constants.ode_params[i]:
                break
        ii = 0
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
        X_ans[-1][1] = 1.0 - X_ans[-1][0] - sum(X_ans[-1][2:])
        X_0 = X_ans[-1]
        i, X_ans, Tres = iter_process(X_0, ts[i], ts, i, X_ans, Tres)
except TechnicalCalcExitException as e:
    i, Tres, X_ans = e.i, e.Tres, e.X_ans
    logging.warning("interrupt by user")


# рисуем всё, что можно :)
time.sleep(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\textbf{t} (s)$')
plt.xlim([constants.to_norm_time(1e-4), constants.to_norm_time(Tres[-1])])
ylabel = r"\textbf{X}"
plt.ylabel(ylabel)
plt.ylim([1e-40, 1000])

new_Tres = np.copy(Tres)
new_Xans = np.copy(X_ans)

for i in range(len(Tres)):
    new_Tres[i] = constants.to_norm_time(new_Tres[i])


register.registrator.calc_plot(plt, new_Tres, X_ans)

for t in new_Tres:
    plt.axvline(x=t, linewidth=0.1)  


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=20)
plt.xscale('log')
plt.yscale('log')

man_aval = {
    "Li_7": (4.885, False),
    "He_3": (0.089, False),
    "Be_7": (0.0964, False),
    "n": (0.005, False)
}

enouth_log_width = 3
not_distinguishable = 1
aval_t = [[0, 0, False] for _ in range(len(X_ans[-1]))]
for i in range(len(Tres) - 20, 0, -1):
    log_x = np.copy(np.log(X_ans[i]))
    log_x_s = sorted(log_x)
    this_up = [0] * len(log_x)
    this_down = [0] * len(log_x)
    for j in range(len(log_x)):
        for k in range(len(log_x_s)):
            if log_x[j] == log_x_s[k]:
                break
        if k == 0:
            this_down[j] = 1000
        else:
            this_down[j] = log_x_s[k] - log_x_s[k-1]
        if k == len(log_x) - 1:
            this_up[j] = 1000
        else:
            this_up[j] = log_x_s[k+1] - log_x_s[k]
    
    for j in range(len(log_x)):
        for man_pos_el in man_aval:
            if man_pos_el in register.registrator.elements[j].names:
                if Tres[i] <= man_aval[man_pos_el][0] < Tres[i+1]:
                    aval_t[j][2] = True
                    aval_t[j][0] = i 
                    aval_t[j][1] = X_ans[i][j] * (1.8 if man_aval[man_pos_el][1] else 0.022)

        if X_ans[i][j] < 1e-40:
            continue
        if aval_t[j][2] == True:
            continue
        if this_down[j] >= not_distinguishable:
            if this_up[j] >= enouth_log_width:
                aval_t[j][2] = True
                aval_t[j][1] = X_ans[i][j]*1.8
                aval_t[j][0] = i
            else:
                if this_up[j] > aval_t[j][1]:
                    aval_t[j][0] = i
                    aval_t[j][1] = X_ans[i][j] * 1.8
            if this_up[j] >= not_distinguishable:
                if this_down[j] >= 2*enouth_log_width:
                    aval_t[j][2] = True
                    aval_t[j][1] = X_ans[i][j] * 0.022
                    aval_t[j][0] = i

        elif this_up[j] >= not_distinguishable:
            if this_down[j] >= 2*enouth_log_width:
                aval_t[j][2] = True
                aval_t[j][1] = X_ans[i][j] * 0.022
                aval_t[j][0] = i

for i in range(len(aval_t)):
    plt.text(new_Tres[aval_t[i][0]], aval_t[i][1], "${}$".format(register.registrator.elements[i].names[-1]))


# plt.title(r"$\eta$= {}".format(constants.nu_n))
plt.xlabel(r'\textbf{time} (s)')
plt.ylabel(r'\textbf{$X_i$}')
# plt.legend(bbox_to_anchor=(1.01, 0.9), loc=2, borderaxespad=0.)
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=2, mode="expand", borderaxespad=0.)
logging.info(("TIME WORKS", (time.time() - start_time)/60))
# plt.title(now_title)
with open("Output.pickle", "wb") as f:
    pickle.dump((X_ans, Tres), f)

sys.stdout.flush()

with open('bbn_with_eta{0:.2f}e-10.csv'.format(constants.nu_n*10**10), 'w', newline='') as csvfile:
    my_len = 7
    mywriter = csv.writer(csvfile, delimiter=';',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    head = ["T9"]+["log10[X({})]".format(el.names[0]) for el in register.registrator.elements]
    for i in range(len(head)):
        head[i] = head[i] + " " * (my_len-len(head[i]))
    mywriter.writerow(head)
    for i in range(len(Tres)):
        try:
            mywriter.writerow(
                ["{0:.8f}".format(constants.to_norm_tempreture(
                    tempreture.Tfromt(Tres[i])
                    , units="T9"))] + 
                ["{0:.8f}".format(math.log10(el)) for el in X_ans[i]]
                )
        except ValueError as e:
            pass

plt.show()
