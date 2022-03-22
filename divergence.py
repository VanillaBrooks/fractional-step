from main import Dims
import numpy as np
from main import NpArray, iu, iv, ip

# divergence maps the velocity to the pressure field
def div(dims: Dims, q: NpArray) -> NpArray:
    Dq = np.zeros((dims.np, 1))

    Dp = np.zeros((dims.np, 1))

    # inner domain

    for i in range(1, dims.nx - 1):
        for j in range(2, dims.ny - 1):
            du_dx = (q[iu[i,j]] - q[iu[i-1, j]]) / dims.dx
            dv_dy = (q[iv[i,j]] - q[iv[i-1, j]]) / dims.dy

            # we are working with velocity here, but the size after the
            # graident should be the same number of points as the pressure field
            # so we index with ip() here
            Dq[ip(i,j)] = du_dx + dv_dy

    # bottom left
    # we need to pin this value
    # we could even pin the divergence at a differnt place like top left or top right
    # avoid top left and top right since it is moving
    Dq[1,1] = 0;

    # bottom right
    # 1 for matlab
    j = 0
    # dims.nx for matlab
    i = dims.nx -1

    du_dx = (u_r[iu[i,j]] - q[iu[i-1, j]]) / dims.dx
    dv_dy = (q[iv[i,j]] - q[iv[i-1, j]]) / dims.dy

    # we are working with velocity here, but the size after the
    # graident should be the same number of points as the pressure field
    # so we index with ip() here
    Dq[ip(i,j)] = du_dx + dv_dy

    # bottom inner

    j = 0
    for i in range(1, dims.nx-1):
        pass

    # also do a bunch of other stuff here


