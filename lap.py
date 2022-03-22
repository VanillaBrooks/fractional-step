# calculates the laplacian
def lap(q: NpArray, iu: NpArray, iv: NpArray, dim: Dims, nq: int) -> NpArray:
    # things for us to do
    #   indexing (iu, iv)
    #       iu comes from iu(i,j)
    #   indicies you are looping over  (i,j)
    #   boundaries
    #       u_r, u_l, u_t, u_b, v_l, v_r, v_b, v_t

    Lq = np.zeros((nq, 1))

    #
    # inner domain
    #

    # matlab: 2:nx-2
    # the indexing here is  nx-2 becuase we dont want the central difference
    # to access a point on the boundary. The y indexing is only to ny-1
    # because even at indexing on the boundary, you are only indexing on the boundary
    # of the x velocity in the y direction, which is not a boundary condition
    for i in range(1, dim.nx - 2):
        # matlab: 2: ny-1
        for j in range(1, dim.ny - 1):
            right = q[iu[i + 1, j]]
            middle = q[iu[i, j]]
            left = q[iu[i - 1, j]]

            top = q[iu[i, j + 1]]
            middle = q[iu[i, j + 1]]
            bottom = q[iu[i, j - 1]]

            x_der = (right - 2 * middle + left) / dim.dx ** 2
            y_der = (top - 2 * middle + bottom) / dim.dy ** 2

            Lq[iu[i, j]] = x_der + y_der

    # matlab: 2:ny-1
    for i in range(1, dim.nx - 1):
        # matlab: 2:ny-2
        for j in range(1, dim.ny - 2):
            right = q[iv[i + 1, j]]
            middle = q[iv[i, j]]
            left = q[iv[i - 1, j]]

            top = q[iv[i, j + 1]]
            middle = q[iv[i, j + 1]]
            bottom = q[iv[i, j - 1]]

            x_der = (right - 2 * middle + left) / dim.dx ** 2
            y_der = (top - 2 * middle + bottom) / dim.dy ** 2

            Lq[iv[i, j]] = x_der + y_der

    ##########
    ########## bottom left corner
    ##########
    i = 0
    j = 0
    # matlab indexing: 1,1
    right = q[iu[i + 1, j]]
    middle = q[iu[i, j]]
    # u_l EXISTS here since its just the boundary
    left = u_l[j]  # this point is exactly on the boundary

    top = q[iu[i, j + 1]]
    middle = q[iu[i, j + 1]]
    # have to create a "ghost node" here because the BC is
    # not perscribed on the exact bottom left corner, so we have to
    # get a value exactly there
    # ---- eq
    # bounary = (i + ghost node) / 2
    # so
    # (2 * boundary - value at i) = ghost node
    # -----eq
    # this is because the node _below_ the bottom left corner
    # does not actually exist on the grid
    bottom = 2 * u_b[i] - q(iu[i, j])

    x_der = (right - 2 * middle + left) / dim.dx ** 2
    y_der = (top - 2 * middle + bottom) / dim.dy ** 2

    Lq[iu[i, j]] = x_der + y_der

    # unsure about this code below here to the next header,
    # it was a really fast copy paste

    right = q[iv[i + 1, j]]
    middle = q[iv[i, j]]
    # have to create a "ghost node" here because the BC is
    # not perscribed on the exact bottom left corner, so we have to
    # get a value exactly there
    left = 2 * u_b[i] - q(iv[i, j])

    top = q[iu[i, j + 1]]
    middle = q[iu[i, j + 1]]
    # this value just exists now, we dont have to worry about it at all
    bottom = v_b[i]

    x_der = (right - 2 * middle + left) / dim.dx ** 2
    y_der = (top - 2 * middle + bottom) / dim.dy ** 2

    Lq[iu[i, j]] = x_der + y_der

    ##########
    ########## bottom right corner
    ##########
    # matlab indcies i = 1, j = nx-1
    # need to use ghost node at the bottom right boundary

    ##########
    ########## bottom inner
    ##########
    j = 1
    # matlab: 2:nx-2
    # blue horizontal triangles as low as you can in the diagram
    # that are not the boundary
    for i in range(1, dim.nx - 2):
        right = q[iu[i + 1, j]]
        middle = q[iu[i, j]]
        left = q[iu[i - 1, j]]

        top = q[iu[i, j + 1]]
        middle = q[iu[i, j + 1]]
        # this value is not actually in the boundaries
        bottom = 2 * u_b[i] - q[iu[i, j]]

        x_der = (right - 2 * middle + left) / dim.dx ** 2
        y_der = (top - 2 * middle + bottom) / dim.dy ** 2

        Lq[iu[i, j]] = x_der + y_der

    # matlab: 2:nx-1
    # red vertical triangles as low as you can on the diagram
    # that are not the boundary
    for i in range(1, dim.nx - 1):
        right = q[iv[i + 1, j]]
        middle = q[iv[i, j]]
        left = q[iv[i - 1, j]]

        top = q[iv[i, j + 1]]
        middle = q[iv[i, j + 1]]
        # this value is just on the boundary
        bottom = v_b[i]

        x_der = (right - 2 * middle + left) / dim.dx ** 2
        y_der = (top - 2 * middle + bottom) / dim.dy ** 2

        Lq[iv[i, j]] = x_der + y_der

    # top left corner

    # top right corner

    # top inner

    # right inner

    # left inner

    return Lq
