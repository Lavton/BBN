"""
Красота и чистота кода немного пренесена в жертву удобству добавления
новых реакций и быстроте считывания.
Достаточно описать реакцию, а ode и jacob построятся сами
"""

import constants
import reactions_list
import logging
from collections import defaultdict, Counter
import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

elist = constants.elist
reactions = reactions_list.reactions
for el in elist:
    exec("import elements.{}".format(el[0]))
    exec("from elements.{0} import {0}".format(el[0]))

def get_num(name):
    for j in range(len(elist)):
        if name in elist[j]:
            return j
    return 0

def get_seq(el):
    num = get_num(el)
    return "(X[{}]/{}.A)".format(num, elist[num][0])

dXss = defaultdict(list)
# для всех реакций
for re in reactions:
    for r in re:
        # "D + He_3 to He_4 + p, elements.He_4.dhe3_he4p, elements.He_4.he4p_dhe3"
        # разбивается на:
        # D + He_3 to He_4 + p
        # elements.He_4.dhe3_he4p
        # elements.He_4.he4p_dhe3
        formula, forw, backw = map(lambda s: s.strip(), r.split(","))
        # D + He_3 to He_4 + p
        # разбивается на:
        # D + He_3
        # He_4 + p
        left, right = map(lambda s: s.strip(), formula.split("to"))
        # D + He_3
        # разбивается на:
        # D
        # He_3
        lefts = list(map(lambda s: s.strip(), left.split("+")))
        # аналогично для He_4 + p
        rights = list(map(lambda s: s.strip(), right.split("+")))
        # считается количество одинаковых элементов для их правильного учёта
        c_left = Counter(lefts)
        c_right = Counter(rights)

        # объединяются элементы, ответственные за распад при прямой реакции
        forw_part = []
        for l in lefts:
            forw_part.append(l)
        forw_part = (forw_part, forw + "(T)/"+str(max(c_left.values())))
        # объединяются элементы, ответственные за синтез при прямой реакции
        back_part = []
        for l in rights:
            back_part.append(l)
        back_part = (back_part, backw + "(T)/"+str(max(c_right.values())))
        
        # для прямой реакции распад идёт со знаком "-"
        # синтез со знаком "+"
        for left in c_left:
            dXss[get_num(left)].append((" - ", c_left[left], forw_part))
            dXss[get_num(left)].append((" + ", c_left[left], back_part))
        # для обратной реакции - распад с "+", синтез с "-"
        for right in c_right:
            dXss[get_num(right)].append((" - ", c_right[right], back_part))
            dXss[get_num(right)].append((" + ", c_right[right], forw_part))

dX = [[] for _ in range(max(dXss.keys()) + 1)]
for k in dXss:
    dX[k] = dXss[k]

def to_ode(seq):
    return "*".join(map(get_seq, seq))

def to_jacob(seq):
    c_seq = Counter(seq)
    res = []
    for c_s in c_seq:
        part = []
        for el in c_seq:
            coef = c_seq[c_s] if el==c_s else 1
            power = coef - 1 if el==c_s else 1
            n_line = "({}*{}**{})".format(coef, get_seq(el), power)
            if el==c_s:
                n_line += "*(1.0/{}.A)".format(elist[get_num(el)][0])
            part.append(n_line)
        res.append((c_s, "*".join(part)))    
    return res

equation = []
for i in range(len(dX)):
    lines = "0"
    for j in range(len(dX[i])):
        lines += "{sign}{coef}*{seq}*{rate}".format(
            sign=dX[i][j][0], 
            coef=dX[i][j][1],
            seq=to_ode(dX[i][j][2][0]),
            rate=dX[i][j][2][1]
        )
    equation.append("my_dx[{}] = ({})*{}.A".format(i, lines, elist[i][0]))
    

ode_line = (
"""
def magic_ode(X, T):
    my_dx = [0 for _ in range(len(X))]
    
    {}
    
    return my_dx
""".format("\n    ".join(equation))
    )

exec(ode_line)

equations = []
for i in range(len(dX)):
    equ_jac = ["0" for _ in range(len(dX))]
    for j in range(len(dX[i])):
        pre_jac = to_jacob(dX[i][j][2][0])
        for prj in pre_jac:
            equ_jac[get_num(prj[0])] += "{sign}{coef}*{seq}*{rate}*{el}.A".format(
            sign=dX[i][j][0], 
            coef=dX[i][j][1],
            seq=prj[1],
            rate=dX[i][j][2][1],
            el=elist[i][0]
        )
    equations.append(equ_jac)

jac_equations = []
for i in range(len(equations)):
    for j in range(len(equations[i])):
        jac_equations.append("my_j[{}][{}] = {}".format(i, j, equations[i][j]))

jacob_line = (
"""
def magic_jacob(X, T):
    my_j = [[0 for _ in range(len(X))] for _ in range(len(X))]
    {}
    
    return my_j
""".format("\n    ".join(jac_equations))
    )

exec(jacob_line)


all_rates = defaultdict(lambda: [[], []])
def set_rates():
    global all_rates
    for re in reactions:
        for r in re:
            # "D + He_3 to He_4 + p, elements.He_4.dhe3_he4p, elements.He_4.he4p_dhe3"
            # разбивается на:
            # D + He_3 to He_4 + p
            # elements.He_4.dhe3_he4p
            # elements.He_4.he4p_dhe3
            formula, forw, backw = map(lambda s: s.strip(), r.split(","))
            # D + He_3 to He_4 + p
            # разбивается на:
            # D + He_3
            # He_4 + p
            left, right = map(lambda s: s.strip(), formula.split("to"))
            # D + He_3
            # разбивается на:
            # D
            # He_3
            lefts = list(map(lambda s: s.strip(), left.split("+")))
            # аналогично для He_4 + p
            rights = list(map(lambda s: s.strip(), right.split("+")))
            # считается количество одинаковых элементов для их правильного учёта
            c_left = Counter(lefts)
            c_right = Counter(rights)
            for right in c_right:
                r = elist[get_num(right)][0]
                all_rates[r][0].append(eval(forw))
                # s1 = "{}.forward_rates.append({})".format(r, forw)
                # exec(s1)
                # print(s1)
                # s2 = "{}.backward_rates.append({})".format(r, backw)
                # exec(s2)
                # print(s2)
                all_rates[r][1].append(eval(backw))
            for left in c_left:
                l = elist[get_num(left)][0]
                all_rates[l][1].append(eval(forw))
                all_rates[l][0].append(eval(backw))
                # exec("{}.forward_rates.append({})".format(l, backw))
                # exec("{}.backward_rates.append({})".format(l, forw))

set_rates()
# print(all_rates)
if __name__ == '__main__':
    print("ODES")
    print(ode_line)
    print()
    print()
    print("JACOB")
    print(jacob_line)