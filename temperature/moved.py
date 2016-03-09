import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import integrate
from scipy.misc import derivative

t_e = 4.3694
m_e = 511000*1.6*10**(-12) # в эргах
k_b = 1.38*10**(-16) # эрг/К

def _s_beaut_integrand_f(x, z):
	if z == 1 or z == 0:
		return 0
			
	try:
		res = (y**2)*(math.sqrt(y**2+x**2)+(y**2)/(3*math.sqrt(y**2+x**2)))*1/(math.exp(math.sqrt(y**2+x**2)+1))
	except OverflowError:
		res = 0.0
	# print (res)
	return res

def S_beaut(x, max_y=np.inf):
	"""3.1.18"""
	integ = scipy.integrate.quad(lambda y: _s_beaut_integrand_f(x, y), 0, max_y)
	return 1 + (45/(2*math.pi**4))*integ[0]


def _eps_beaut_integrand_f(x, y):
	try:
		res = y**2*math.sqrt(y**2+x**2)/(math.exp(y**2+x**2)+1)
	except OverflowError:
		res = 0.0
	return res

def epsilon(x, max_y=np.inf):
	"""3.1.23"""
	integ = scipy.integrate.quad(lambda y: _eps_beaut_integrand_f(x, y), 0, max_y)
	return 1+(21/8.)*(4.*S_beaut(x)/11)**(4./3)+(30/(math.pi**4))*integ[0]

def under_int(x):
	return (3-x*derivative(S_beaut, x)/S_beaut(x))*(epsilon(x)**(-1/2.))*x

def tfromT(T):
	"""3.1.24"""
	integ = t_e*scipy.integrate.quad(under_int, 0, m_e/(k_b*T))[0]
	return integ


if __name__ == '__main__':
	x = 1.0
	t0 = tfromT(10**11)
	Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=100)
	# y = np.logspace(math.log(0.1), math.log(10.))
	# xs = np.logspace(math.log(0.1), math.log(100.), num=100)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'\textbf{tempreture} (K)')
	plt.ylabel(r'\textbf{time} (s)')
	plt.plot(Ts, [0.994*(T/10**10)**(-2)-t0 for T in Ts], 
		linewidth=2.0, label=r't\to 0')
	plt.plot(Ts, [1.78*(T/10**10)**(-2)-t0 for T in Ts], 
		linewidth=2.0, label=r't\to $\infty$')
	plt.plot(Ts, [tfromT(T)-t0 for T in Ts], 
		'r--', label=r't(T) modeling result')
	plt.legend()
	plt.show()
	print([0.994*(T/10**10)**(-2) for T in Ts])
	# print(y)
	# y = [1., 2., 10., 100., 500.]
	# plt.plot(y, [_s_beaut_integrand_f(x,ys) for ys in y])
	# print([S_beaut(x,ys) for ys in y])
	# plt.plot(y, [S_beaut(x,ys) for ys in y])
	# plt.show()
	# print([S_beaut(x) for x in xs])
	# print([epsilon(x) for x in xs])
	# plt.plot(xs, [S_beaut(x) for x in xs])
	# plt.show()
	# plt.plot(xs, [epsilon(x) for x in xs])
	# plt.show()
	T = 10**10
	print("{:.1E}: {:.4E}".format(T, tfromT(T)-t0))