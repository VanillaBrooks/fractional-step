from main import Dims
import numpy as np
from main import NpArray, iu, iv

def adv(q: NpArray, dim: Dims) -> NpAray:

    ##########
    ########## Inner domain
    ########## 

    Np = np.zeros((dims.nq, 1))
    # Nq = diff(u_x  * u_x, x) + diff (u_y * v_x, y);

    # matlab: 2:ny-2
    for i in range(1,dim.nx-2):
        # matlab: 2:ny-1
        for j in range(1,dim.ny-1):
            # this is just a simple forward difference
            # doing a central difference makes this life a bit more complex
            u_x_i_plus_half = (1/2) * (q[iu[i+1, j]] + q[iu[i,j]])
            u_x_i_minus_half = (1/2) * (q[iu[i, j]] + q[iu[i-1,j]])

            diff_u_x = (u_x_i_plus_half - u_x_i_minus_half)/dims.dx

            # this might be off with indexing / iv stuff
            u_y_plus_half = (1/2) * (q[iu[i, j+1]] + q[iu[i,j]])
            u_y_minus_half = (1/2) * (q[iu[i, j]] + q[iu[i,j-1]])

            v_x_plus_half= (1/2) * (q[iv[i+1, j]] + q[iv[i,j]])
            v_x_minus_half= (1/2) * (q[iv[i, j]] + q[iv[i-1,j]])

            diff_second_term = (u_y_plus_half * v_x_plus_half - u_y_minus_half * u_y_plus_half)/dims.dy

            # missing stuff

            Nq(iu[i,j]) = (term**2 + forward**2) / dims.dx

    
    # matlab: 2:ny-1
    for i in range(1,dim.nx-1):
        # matlab: 2:ny-2
        for j in range(1,dim.ny-2):
