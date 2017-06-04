import numpy
import logging
"""fundamental constants"""

ergToEV = 6.2415 * 10**11 # перевод из эргов в электроны вольты
eVToErg = 1.6 * 10**(-12) # перевод из электрон-вольтов в эрги
m_e = 511000 # в эВ
k_b = 1.38 * 10**(-16) * ergToEV # эВ/К
Q = 1.293 * 10**6 # разница между протоном и нейтроном, эВ
initT = 10**11 # начальная температура, К
g = 1.6638 * 10**4 # константа взаимодействия Ферми, эВ^-2
c = 29979245800 # скорость света, cм/с
h = 6.582 * 10**(-16) # приведённая постоянная планка, эВ/с
t_n = 880 # в с, время жизни нейтрона
nu_n = 6.1 * 10**-10 # барион-фотонное отношение
nu_0 = 6 * 10**-10
m_p = 938.272 * 10**6 # масса протона, эВ
m_n = 939.565 * 10**6 # масса нейтрона, эВ
N_a = 6 * 10**23 # число авагадро
m_d = 2.01355321271 * 1.66053873e-24*29979245800.0*29979245800.0 # масса дейтерия, эВ
amuToErg = 1.66053873e-24*c*c

########################################################
# технические константы
sql_enabled = True # разрешить ли кеширование в БД
sql_enabled = sql_enabled if sql_enabled else False
smart_caching = True # запускает продолжение кода с момента остановки/изменеия параметров


func_name_to_db_name = {
    "tfromT": "TimeFromTempreture",
    "Tfromt": "TempretureFromTime",
    "TnuFromT": "TempreNuFromTempreture",
    "lambda_n__p": "LambdaNPFromTempr",
    "lambda_p__n": "LambdaPNFromTempr",
    "derriviate_T_from_t": "DerriviateTdtFromTempreture"
}


ode_params = [
    [-1.0, {
    "rtoi": 1e-6,
    "max_step": 1.0,
    "min_step": 1e-10,
    }], 
    [0.00045, {
    "atoi": 1e-10,
    "rtoi": 1e-10,
    "min_step": 0.0,
    }
    ],
    [1.29, {
        "max_step": 0.01,
        "rtoi": 1e-9,
        "atoi": 1e-9,
    }
    ],
    [300, {
    "max_step": 1.0,
    },
    [1000, {
        "max_step": 0.0
    }]
    ]
]

def_params = {
    "rtoi": 1e-6,
    "atoi": 1e-12,
    "max_step": 1.0,
    "min_step": 1e-10,
    "last_step": 0.0
}

prog_status = "DEBUG"


elist = [
    ["n"],
    ["H_1", "p"],
    ["H_2", "d", "D"],
    ["He_3"],
    ["H_3", "t", "T"],
    ["He_4"],
    ["Be_7"],
    ["Li_7"],
    ["Li_6"],
    ["He_6"]
]



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

    T_ *= m_n
    if units=="MeV":
        T_ /= 10**6
    if units=="K":
        T_ = T_ / k_b
    if units=="T9":
        T_ = (T_ / k_b) * 10**(-9)
    return T_

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

E_d_t0=less_tempreture(m_p+m_n-m_d)
