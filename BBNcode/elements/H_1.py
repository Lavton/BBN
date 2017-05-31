"""протон, или $^{1}H$"""

import constants
import math
import nTOp
import functools
from elements.Element import Element

n_to_p = nTOp.lambda_n__p
p_to_n = nTOp.lambda_p__n

# создать новый элемент
H_1 = Element("H_1", 0.5)

H_1.forward_rates.append(nTOp.lambda_n__p)
H_1.backward_rates.append(nTOp.lambda_p__n)