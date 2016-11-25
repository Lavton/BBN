import numpy
"""fundamental constants"""

ergToEV = 6.2415 * 10**11 # перевод из эргов в электроны вольны
eVToErg = 1.6 * 10**(-12) # перевод из электрон-вольтов в эрги
m_e = 511000 # в эВ
k_b = 1.38 * 10**(-16) * ergToEV # эВ/К
Q = 1.293 * 10**6 # разница между протоном и нейтроном, эВ
initT = 10**11 # начальная температура, К
g = 1.6638 * 10**4 # константа взаимодействия Ферми, эВ^-2
c = 29979245800 # скорость света, cм/с
h = 6.582 * 10**(-16) # приведённая постоянная планка, эВ/с
t_n = 880 # в с, время жизни нейтрона
nu_n = 6 * 10**-10 # барион-фотонное отношение
m_p = 938 * 10**6 # масса протона, эВ
m_n = 939.565 * 10**6 # масса нейтрона, эВ
N_a = 6 * 10**23 # число авагадро



########################################################
# технические константы
sql_enabled = True # разрешить ли кеширование в БД
sql_enabled = sql_enabled if sql_enabled else False

func_name_to_db_name = {
    "tfromT": "TimeFromTempreture",
    "TnuFromT": "TempreNuFromTempreture",
    "lambda_n__p": "LambdaNPFromTempr",
    "lambda_p__n": "LambdaPNFromTempr"
}

#########################################################
##### преобразования в размерную величину ###############

def less_tempreture(T, units="eV"):
    T_ = T
    if type(T) is numpy.ndarray:
        T_ = T.copy()
    if units=="K":
        T_ = k_b * T_
    return T_ / m_n

def to_norm_tempreture(T, units="eV"):
    T_ = T
    if type(T) is numpy.ndarray:
        T_ = T.copy()
    if units=="K":
        T_ = T_ / k_b
    return T_ * m_n

def less_time(t):
    t_ = t
    if type(t) is numpy.ndarray:
        t_ = t.copy()
    return t_ / t_n

def to_norm_time(t):
    t_ = t
    if type(t) is numpy.ndarray:
        t_ = t.copy()
    return t_ * t_n