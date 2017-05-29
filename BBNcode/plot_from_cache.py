import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pickle

if not os.path.isfile("smart_cache.pickle"):
    exit()


with open("smart_cache.pickle", "rb") as f:
    i, X_ans, Tres, ode_params = pickle.load(f)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\textbf{t} (s)$')
plt.xlim([1e-4,1e3])
ylabel = r"\textbf{X}"
plt.ylabel(ylabel)
plt.ylim([1e-30, 1000])
import elements.register as elements
elements.registrator.calc_plot(plt, Tres, X_ans)

for t in Tres:
  plt.axvline(x=t, linewidth=0.1)

from elements.H_2 import H_2
plt.axvline(x=H_2.tr_t)
from elements.He_3 import He_3
plt.axvline(x=He_3.tr_t)
from elements.H_3 import H_3
plt.axvline(x=H_3.tr_t)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'\textbf{time} (s)')
plt.ylabel(r'\textbf{\lambda}')
plt.legend()
plt.show()
