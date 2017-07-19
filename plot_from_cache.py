import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pickle
import numpy as np 
import tempreture

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params, element_stuct = pickle.load(f)

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
    rels = []
    for el in element_stuct:
        if el[1]:
            rels.append((el[1], "${}$".format(el[0][-1])))
    rels.sort()
    start = 1e-4
    step = 1e5
    for el in rels:
      plt.axvline(x=el[0])
      plt.text(el[0], start, el[1])  
      start /= step
    for op in ode_params:
        plt.axvline(x=op[0], linewidth=0.5, color="red")

print_addition()

import register
register.registrator.calc_plot(plt, Tres, X_ans)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{time} (s)')
plt.ylabel(r'\textbf{\lambda}')
plt.legend(bbox_to_anchor=(1.01, 0.9), loc=2, borderaxespad=0.)
plt.show()
