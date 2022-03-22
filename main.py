import numpy as np
from adv import adv
from lap import lap

NpArray = np.ndarray

# QUESTIONS
# 1) are ip, iv, iu matricies that map 2d coordinates to their 1D counterparts? 
#   I thought they were functions


class Dims:
    def __init__(self, cavity_length: float, cavity_height: float):
        # number of cells in the x direction
        self.nx = 32
        # number of cells in the y direction
        self.ny = 32

        # number of points of pressure
        # (since pressure is
        self.np = self.nx * self.ny

        # velocity points (say in the x direction) have all
        # the points in the y direction, but only nx-1
        # points in the x direction (since we dont explicity
        # handle values at the bounary)
        # we include both of these for the vector lengths
        # since the velocity vector will include both
        self.nu = (self.nx - 1) * self.ny + (self.ny - 1) * self.nx

        # grid spacing in x and y directions
        self.dx = cavity_length / self.nx
        self.dy = cavity_height / self.ny


# Take an array and the 2D indexing of the array
# and return the linear index at which that value is
def to_linear_index(i: int, j: int) -> int:
    pass

def main():

    cavity_length = 1
    cavity_height = 1

    # init the dimensions and
    dim = Dims(cavity_length, cavity_height)

    # Reynolds number
    re = 1000.0

    # number of timesteps
    n_step = 1000

    # vector values of the quantities we are solving for
    # Note from later: what the fuck 
    iu = np.zeros((dim.nx - 1, dim.ny))
    iv = np.zeros((dim.nx, dim.ny - 1))
    ip = np.zeros((dim.nx, dim.ny))

    # velocity field vectors
    q_star = np.zeros((dim.nu, 1))

    # boundary conditions
    # EXERCISE: What are the dimensions?

    # top
    u_t = np.ones((1, 1))  # ones since this is lid drivien cavity
    v_t = np.zeros()

    # bottom
    u_b = np.zeros()
    v_b = np.zeros()

    # left
    u_l = np.zeros()
    v_l = np.zeros()

    # right
    u_r = np.zeros()
    v_r = np.zeros()

    #
    # Assign time stepping variables
    #
    umax = np.max(u_t)  # comes from moving the boundary at the initial velocity
    # there is something called the
    # and another one called
    cfl = 1

    # advective cfl  = U dt / dx
    dt_advective = cfl * min(dim.nx, dim.ny) / umax

    # diffusion cfl = dt / (Re * dx^2)
    dt_diffusive = re * min(dim.nx, dim.ny) ** 2 * cfl

    dt = min(dt_diffusive, dt_advective)

    #
    # variables for fractional step
    #

    # ----- velocity data -----

    # previous (n-1) data
    q_nm1 = q_star
    # current (n) data
    q_n = q_star
    # future (n+1) data
    q_np1 = q_star

    # ----- pressure data -----

    # current (n)
    p_n = np.zeros((dim.np, 1))
    # future (n+1)
    p_np1 = np.zeros((dim.np, 1))

    for iter in range(0, n_step):
        q_nm1 = q_n
        q_n = q_np1

        p_n = p_np1

        # iu(:, :) is a function that takes in two points and outputs their linear indexing for x velocity
        # iv(:, :) is the same, but for V velocity

        Sq_rhs = (
            1 / dt * q_n
            + 1 / (2 * Re) * lap(q_n)
            + (3 / 2) * adv(q_n)
            - (1 / 2) * adv(q_nm1)
        )
        Rq_lhs = 1 / dt * q_n - 1 / (2 * Re) * lap_nbc(q_n)
        # for lap_nbc do the same as lap, but just set all the boundary values
        # to zero, including ghost node values
        q_star = conjugate_gradient(Rq_lhs, Sq_rhs, tol)

        #
        # second fractional step
        #
        lhs = pass
        rhs = pass
        p_np1 = conjugate_gradient(lhs, rhs, tol)

        # third fractional step

        #
        # TODO: code up the three steps of the fractional step method
        #

        #
        # Goals for thursday
        #   - linear indexing operator
        #   - boundarey conditions
        #   - conjugate gradient method

if __name__ == "__main__":
    main()
