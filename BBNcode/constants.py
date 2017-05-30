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
nu_n = 6 * 10**-10 # барион-фотонное отношение
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
    # [0.00046, {
    #     "min_step": 0.0,
    #     "atoi": 1e-8,
    #     "rtoi": 1e-8,
    #     "max_step": 1e-7,
    # }],
    # [0.00003, {
    #     "min_step": 0.0,
    #     "atoi": 1e-1,
    #     "rtoi": 1e-1,
    #     # "max_step": 1e-9,
    # }],
    [0.0013, {
    "atoi": 1e-8,
    "rtoi": 1e-8,
    }
    ],
    [0.0016, {
    "min_step": 0.0,
    "rtoi": 1e-10,
    "atoi": 1e-10,
    }
    ],
    # [0.004, {
    # "rtoi": 1e-7,
    # }],
    [0.09, {
    "rtoi": 1e-11,
    "max_step": 0.0005
    ,
    }],
    [0.02, {
    "atoi": 1e-9,
    "rtoi": 1e-11,
    "max_step": 0.0001,
    "min_step": 0.0
    }],
    [0.038, {
    "atoi": 1e-9,
    "rtoi": 1e-8,
    "max_step": 0.002,
    }],
    [0.052, {
    "rtoi": 1e-6,
    "max_step": 0.2,
    "min_step": 0.0
    }],
    [0.055,{
    "rtoi": 1e-9,
    "max_step": 0.2,
    "atoi": 1e-8,
    }],
    [0.080,{
    "rtoi": 1e-7,
    "atoi": 1e-8,
    # "min_step": 1e-8,
    "max_step": 0.03,
    }],
    [0.2, {
    "rtoi": 1e-8,
    "atoi": 1e-10,
    "max_step": 0.02
    # "not_use_jacob": True
    }],
    [0.31, {
    "rtoi": 1e-11,
    "atoi": 1e-12,
    "max_step": 0.3,
    # "not_use_jacob": True
    }],
    [10, {
    "rtoi": 1e-9,
    "max_step": 0.0
    }],
    [100, {
    "max_step": 0.0,
    }]
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
    ["Li_7"]
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
