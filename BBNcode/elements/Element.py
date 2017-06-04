import constants
import math
import nTOp
import functools

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