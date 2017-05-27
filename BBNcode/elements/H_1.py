"""протон, или $^{1}H$"""

import constants
import math
import nTOp
import functools
from elements.Element import Element

n_to_p = nTOp.lambda_n__p
p_to_n = nTOp.lambda_p__n


H_1 = Element("H_1", 0.5)
H_1.ode_elem = {
    "n": (lambda X, T: nTOp.lambda_p__n(T) * X[1]),
    "H_1": (lambda X, T: nTOp.lambda_n__p(T) * X[0] - nTOp.lambda_p__n(T) * X[1])
}

# внешний - столбец, внутренний - строка

H_1.jacob = {
    "n": {
        "H_1": 
        (lambda X, T: +nTOp.lambda_p__n(T))
    },
    "H_1": {
        "n": 
        (lambda X, T: +nTOp.lambda_n__p(T)),
        "H_1": (lambda X, T: -nTOp.lambda_p__n(T))
    }
}

H_1.forward_rates.append(nTOp.lambda_n__p)
H_1.backward_rates.append(nTOp.lambda_p__n)