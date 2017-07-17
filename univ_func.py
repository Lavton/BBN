from constants import *
import constants
import math
from math import pi

def __photon_density__(T):
    """
    плотность фотонов во вселенной
    """
    # T in eV
    zeta_3 = 1.202
    return 0.244*(constants.to_norm_tempreture(T)/(h*c))**3

def __proton_mass_density__(T):
    """плотность протонов во вселенной"""
    return (m_p/5.60958835719e+32) * __photon_density__(T) * constants.nu_n

def rat_scale(T):
    return __proton_mass_density__(T)  

def universe_speed(T):
    return ((8*pi/3)*constants.g_grav * __proton_mass_density__(T))**(1./2)

if __name__ == '__main__':
	print(h, hbar)
	print(__photon_density__(constants.less_tempreture(2.72548, units="K")))