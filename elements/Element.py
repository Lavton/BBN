import constants
import math
import nTOp
import functools
import univ_func

class Element():
    """
    Общий класс для всех элементов.
    у каждого есть начальное значение, строковое представление, массовое число и функции самого подсчёта
    """

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        self.A = 1
        self.equilibrium = None
        self.forward_rates = []
        self.backward_rates = []
        self.tr_T = None
        self.tr_t = None
        self.interpo = {}
        self.reactions = []
        self.names = [str_view]

    def add_interpo(self, f_name, string_to_parse):
        lines = string_to_parse.strip().replace("−", "-").split("\n")
        rates = []
        for i in range(len(lines)):
            lines[i] = lines[i].split()
            rates.append((float(lines[i][0]), float(lines[i][1])))
            if len(lines[i]) > 6:
                rates.append((float(lines[i][4+0]), float(lines[i][4+1])))
        rates.sort()
        from scipy.interpolate import interp1d
        import numpy as np
        xs = []
        ys = []
        for r in rates:
            T_ = constants.less_tempreture(r[0]*10**9, units="K")
            xs.append(math.log(T_))
            ys.append(math.log(r[1]))
        xs = np.array(xs)
        ys = np.array(ys)
        cs = interp1d(xs, ys, kind='cubic')
        self.interpo[f_name] = cs

    def get_inter_val(self, f_name, T):
        return math.exp(self.interpo[f_name]([math.log(T)])[0])

    def show_rates(self):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import numpy as np
        from tempreture import tfromT
        grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
        grid2 = grid
        Ts = constants.less_tempreture(grid2, units="K")
        ts = np.array([tfromT(T) for T in Ts])
        plt.cla()
        plt.clf()
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'\textbf{tempreture} (MeV)')
        plt.ylabel(r'\textbf{\lambda}')
        i = 1
        for f in self.forward_rates:
            plt.plot(ts, [f(T) for T in Ts], label='forw {}'.format(i))
            i += 1
        i = 1
        for b in self.backward_rates:
            plt.plot(ts, [b(T) for T in Ts], label='back {}'.format(i))
            i += 1
        plt.legend()
        plt.show()
        plt.cla()
        plt.clf()

    def equilib_zeroize(self, func):
        """
        Преобразование в 0 до tr
        """
        @functools.wraps(func)
        def inner(*args, **kwargs):
            return func(*args, **kwargs) if args[0] < self.tr_T else 0
        return inner

    def nacreII(self, func):
        """spline"""
        @functools.wraps(func)
        def inner(*args, **kwargs):
            T9 = constants.to_norm_tempreture(args[0], units="T9")
            ro_b = univ_func.rat_scale(args[0])
            return func(*args, **kwargs) if T9 > 10 else self.get_inter_val(func.__name__, args[0]) * ro_b/(constants.less_time(1))
        return inner

    def set_mass_excess(self, total_mass, n_N, p_N):
        m_ns = constants.m_n * n_N
        m_ps = constants.m_p * p_N
        m = total_mass * (10**(-6)) * constants.amuToErg * constants.ergToEV 
        self.mass = m
        self.mass_excess = constants.less_tempreture(-(m - m_ns - m_ps - constants.m_e*p_N), units="eV")

if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)