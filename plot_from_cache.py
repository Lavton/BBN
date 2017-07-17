import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pickle
import numpy as np 
import tempreture

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params = pickle.load(f)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\textbf{t} (s)$')
plt.xlim([1e-5,1e3])
# plt.xlim([0.0009, 0.0017])
ylabel = r"\textbf{X}"
plt.ylabel(ylabel)
plt.ylim([1e-100, 1000])
# plt.ylim([5e-66, 3e-58])

for t in Tres:
    plt.axvline(x=t, linewidth=0.1)

def print_addition():
    from elements.H_2 import H_2
    plt.axvline(x=H_2.tr_t)
    plt.text(H_2.tr_t, 2e-10, "$H_2$")
    from elements.He_3 import He_3
    plt.axvline(x=He_3.tr_t)
    plt.text(He_3.tr_t, 2e-7, "$He_3$")
    from elements.H_3 import H_3
    plt.axvline(x=H_3.tr_t)
    plt.text(H_3.tr_t, 2e-15, "$H_3$")
    from elements.He_4 import He_4
    plt.axvline(x=He_4.tr_t)
    plt.text(He_4.tr_t, 2e-4, "$He_4$")
    from elements.Be_7 import Be_7
    plt.axvline(x=Be_7.tr_t)
    plt.text(Be_7.tr_t, 2e-4, "$Be_7$")
    from elements.Li_7 import Li_7
    plt.axvline(x=Li_7.tr_t)
    plt.text(Li_7.tr_t, 2e-7, "$Li_7$")
    from elements.Li_6 import Li_6
    plt.axvline(x=Li_6.tr_t)
    plt.text(Li_6.tr_t, 5e-6, "$Li_6$")
    from elements.He_6 import He_6
    plt.axvline(x=He_6.tr_t)
    plt.text(He_6.tr_t, 5e-30, "$He_6$")
    for op in ode_params:
        plt.axvline(x=op[0], linewidth=0.5, color="red")
    # r_he = []
    # for i in range(len(Tres)):
    #     print(i, len(Tres))
    #     he_x = np.copy(X_ans[i])
    #     try:
    #         pass
    #         he_x = He_6.equilibrium([he_x], tempreture.Tfromt(Tres[i]))[0]
    #         r_he.append(np.copy(he_x)[9])
    #     except ZeroDivisionError as e:
    #         print("e")
    #         r_he.append(0.5)
    # plt.plot(Tres, r_he, label="eq")

print_addition()

import elements.register as elements
elements.registrator.calc_plot(plt, Tres, X_ans)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{time} (s)')
plt.ylabel(r'\textbf{\lambda}')
plt.legend(bbox_to_anchor=(1.01, 0.9), loc=2, borderaxespad=0.)
plt.show()