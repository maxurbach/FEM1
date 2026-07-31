# Discretization of the square plate

# IDEA: 
# The plate is subdivided into (nx * ny) quadrilaterals
# and each quadrilateral is split into two triangules by a diagonal

# VISUALIZATION: 
# UL = Upper Left, UR = Upper Right, LL = Lower Left, LR = Lower Right
#
# UL -------- UR
# |         / |
# |      /    |
# |    /      |
# | /         |
# LL -------- LR
#
# We always cut from LL to UR

# NODE NUMBERING:
# Node (i, j) gets global number n = j * (nx + 1) + i + 1
# so we sweep in x-direction first (i), then move up in y-direction (j)

# DEGREES OF FREEDOM:
# Every node has 2 DOFs

# RETURNS:
# coords [nnodes * 2] nodal coordinates [x, y]
# conn   [nelem * 3]  node numbers per element
# eft    [nelem * 6]  global DOF numbers per element
# ndof   total number of DOFs

function create_mesh(nx, ny, L = 1.0)
    nnodes = (nx + 1) * (ny + 1)    # number of nodes
    nelem = 2 * nx * ny             # number of traingles (2 per quadrilateral)
    ndof = 2 * nnodes               # 2 DOFs per node (ux, uy)

    # nodal coordinates
    # coords[n, 1] = x, coords[n, 2] = y
    coords = zeros(nnodes, 2)
    for j in 0:ny
        for i in 0:nx
            n = j * (nx + 1) + i + 1    # global node number
            coords[n, 1] = i * L / nx   # x = i * dx
            coords[n, 2] = j * L / ny   # y = j * dy
        end
    end

    # connectivity and EFT
    conn = zeros(Int, nelem, 3)     # node numbers [n1, n2, n3]
    eft = zeros(Int, nelem, 6)      # global DOF numbers [ux1, uy1, ux2, uy2, ux3, uy3]

    e = 0   # element counter
    for j in 0:ny-1
        for i in 0:nx-1
            # node number at the four corners of the current quadrilateral
            nLL = j       * (nx + 1) + i       + 1     # lower left
            nLR = j       * (nx + 1) + (i + 1) + 1     # lower right
            nUL = (j + 1) * (nx + 1) + i       + 1     # upper left
            nUR = (j + 1) * (nx + 1) + (i + 1) + 1     # upper right  

            # triangle 1: LL, LR, UR (lower triangle)
            e += 1
            conn[e, :] = [nLL, nLR, nUR]
            eft[e, :] = [2 * nLL - 1, 2 * nLL, 2 * nLR - 1, 2 * nLR, 2 * nUR - 1, 2 * nUR]

            # triangle 2: LL, UR, UL (upper triangle)
            e += 1
            conn[e, :] = [nLL, nUR, nUL]
            eft[e, :] = [2 * nLL - 1, 2 * nLL, 2 * nUR - 1, 2 * nUR, 2 * nUL - 1, 2 * nUL]
        end
    end

    return coords, conn, eft, ndof
end
