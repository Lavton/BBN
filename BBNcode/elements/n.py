"""
нейтроны
"""

import constants
import math
import nTOp
import functools
from elements.Element import Element

n_to_p = nTOp.lambda_n__p
p_to_n = nTOp.lambda_p__n

n = Element("n", 0.5)
n.ode_elem = {
    "n" : (lambda X, T: -nTOp.lambda_n__p(T) * X[0])
}

n.jacob = {
    "n": {
        "n": (lambda X, T: -nTOp.lambda_n__p(T))
    }
}

n.set_mass_excess(1008664.9157, n_N=1, p_N=0)
n.backward_rates.append(nTOp.lambda_n__p)
n.forward_rates.append(nTOp.lambda_p__n)