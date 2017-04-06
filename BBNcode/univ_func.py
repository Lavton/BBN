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
    # return (2 * (8*pi)/c**3 * zeta_3) / (h/(2*pi*constants.to_norm_tempreture(T)))**3
    return 0.244*(constants.to_norm_tempreture(T)/(h*c))**3

def __proton_mass_density__(T):
    """плотность протонов во вселенной"""
    return (m_p/5.60958835719e+32) * __photon_density__(T) * nu_n

def rat_scale(T):
    return __proton_mass_density__(T)   # * constants.N_a
    # return __photon_density__(T) * nu_n * constants.N_a