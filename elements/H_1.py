"""протон, или $^{1}H$"""

import constants
import math
import nTOp
import functools
from elements.Element import Element

def n_to_p(T):
    return nTOp.lambda_n__p(T)
def p_to_n(T):
    return nTOp.lambda_p__n(T)

# создать новый элемент
H_1 = Element("H_1", 0.5)
H_1.set_mass_excess(1007825.03207, n_N=0, p_N=1)
H_1.forward_rates.append(nTOp.lambda_n__p)
H_1.backward_rates.append(nTOp.lambda_p__n)
H_1.names = ["H_1", "p", "^1H"]
H_1.reactions.append((
    ("p",), 
    ("n",), 
    p_to_n, 
    n_to_p
    ))
    