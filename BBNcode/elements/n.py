"""
нейтроны
"""

import constants
import math
import nTOp

from elements.Element import Element

n = Element("n", 0.5)
n.ode_elem = {
    "n" : (lambda X, T: -nTOp.lambda_n__p(T) * X[0])
}

n.jacob = {
    "n": {
        "n": (lambda X, T: -nTOp.lambda_n__p(T))
    }
}

