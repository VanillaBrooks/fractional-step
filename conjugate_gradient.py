from main import Dims
import numpy as np
from main import NpArray, iu, iv, ip



# if flag = 1
#   we are solving the first step of the fractional step
# if flag = 2
#   solve (second step ?)
# PARAMS:
# lhs -> Ax
# rhs -> b
#
# residual = (Ax - b)
# d = residual' * residual
# 
# A is always accompanied by and x, and they are unique to each other
def conjugate_gradient(lhs: NpArray, rhs: NpArray, tolerance: float, flag: int) -> NpArray:
    if flag == 1:
        # ????
        Rq_lhs = 1 / dt * q_n - 1 / (2 * Re) * lap_nbc(q_n)
    elif flag == 2:
        pass
