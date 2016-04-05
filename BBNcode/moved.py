import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sys
from scipy import integrate
from scipy.misc import derivative

t_e = 4.3694
m_e = 511000*1.6*10**(-12) # в эргах
k_b = 1.38*10**(-16) # эрг/К

def _s_beaut_integrand_f(x, z):
	if z == 1 or z == 0:
		return 0
	y = 1/z - 1
	try:
		res = (y**2)*(math.sqrt(y**2+x**2)+(y**2)/(3*math.sqrt(y**2+x**2)))*1/(math.exp(math.sqrt(y**2+x**2))+1)/z**2
	except OverflowError:
		# print ("ERR: {}".format(y))
		res = 0.0
	return res

def S_beaut(x):
	"""3.1.18"""
	integ = scipy.integrate.quad(lambda z: _s_beaut_integrand_f(x, z), 0, 1)
	return 1 + (45/(2*math.pi**4))*integ[0]

def _eps_beaut_integrand_f(x, z):
	if z == 1 or z == 0:
		return 0
	y = 1/z - 1
	try:
		res = ((y/z)**2)*math.sqrt(y**2+x**2)/(math.exp(math.sqrt(y**2+x**2))+1)
	except OverflowError:
		res = 0.0
	return res

def epsilon(x):
	"""3.1.23"""
	integ = scipy.integrate.quad(lambda z: _eps_beaut_integrand_f(x, z), 0, 1)
	return 1+(21/8.)*(4.*S_beaut(x)/11)**(4./3)+(30/(math.pi**4))*integ[0]

def under_int(x):
	return (3-x*derivative(S_beaut, x)/S_beaut(x))*(epsilon(x)**(-1/2.))*x

def tfromT(T):
	"""3.1.24"""
	integ = t_e*scipy.integrate.quad(under_int, 0, m_e/(k_b*T))[0]
	return integ


if __name__ == '__main__':
	if len(sys.argv)-1:
		T = sys.argv[1]
		T = eval(T)
		print("{:.1E}: {:.4E}\n".format(T, tfromT(T)))
		exit()
	t0 = tfromT(10**11)
	Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=100)
	ts = [tfromT(T)-t0 for T in Ts]
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'\textbf{tempreture} (K)')
	plt.ylabel(r'\textbf{time} (s)')
	plt.plot(Ts, [0.994*(T/10**10)**(-2)-0.994*(10)**(-2) for T in Ts], 
		linewidth=2.0, label=r't\to 0')
	plt.plot(Ts, [1.78*(T/10**10)**(-2)-1.78*(10)**(-2) for T in Ts], 
		linewidth=2.0, label=r't\to $\infty$')
	plt.plot(Ts, ts, 
		'r--', label=r't(T) modeling result')
	plt.legend()
	plt.show()
	with open("tempreture_data.dat", "w") as f:
		for i in range(len(Ts)):
			f.write("{:.1E}: {:.4E}\n".format(Ts[i], ts[i]))