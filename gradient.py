from main import Dims
import numpy as np
from main import NpArray, iu, iv, ip

# gradient maps the pressure to velocity variables
def gradient(dims: Dims) -> NpArray:
    # there are no boundaries to compute
    # because R is calculated with the laplacian L that
    # does not calculate values at boundaries
    
    grad = np.zeros((dims.nu, 1))

    # dp/dx comp

    #2 - nx for matlab
    for i in range(1, dims.nx):
        for j in range(1, dims.ny):
            grad[iu(i,j)] = (p[ip[i,j]] - p[ip[i-1, j]]) / dims.dx

    # dp/dy component
    # 1: nx
    for i in range(0, dims.nx):
        # 2: ny
        for i in range(1, dims.ny):
            pass
            # something else here
