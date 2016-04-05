import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sys
from scipy import integrate
from scipy.misc import derivative
from tempreture import *


if __name__ == '__main__':
    if len(sys.argv)-1:
        T = sys.argv[1]
        T = eval(T)
        print("{:.1E}: {:.4E}\n".format(T, tfromT(T)))
        exit()
    Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=100)
    ts = [tfromT(T) for T in Ts]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{time} (s)')
    plt.ylabel(r'\textbf{tempreture} (K)')
    plt.plot([0.994*(T/10**10)**(-2)-0.994*(10)**(-2) for T in Ts], Ts, 
        linewidth=2.0, label=r't\to 0')
    plt.plot([1.78*(T/10**10)**(-2)-1.78*(10)**(-2) for T in Ts], Ts,
        linewidth=2.0, label=r't\to $\infty$')
    plt.plot(ts, Ts,
        'r--', label=r't(T) modeling result')
    plt.legend()
    plt.show()
    with open("tempreture_data.dat", "w") as f:
        for i in range(len(Ts)):
            f.write("{:.1E}: {:.4E}\n".format(Ts[i], ts[i]))