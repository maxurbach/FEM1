# DESCRIPTION:
# 
# Generates a structured finite element mesh of 4-node bilinear isoparametric quadrilateral 
# elements (Quad4) for a 90-degree curved beam (quarter annulus) using polar coordinate 
# mapping (r, theta) -> (x, y)

# INPUT:
# nr     - number of elements in the radial direction
# ntheta - number of elements in the circumferential direction
# ri     - inner radius of the beam
# ro     - outer radius of the beam

# OUTPUT:
# coords - [nnodes * 2] nodal coordinates [x1, x2] <- coordinate system from figure 1
# conn   - [nelem * 4] node numbers per element
# eft    - [nelem * 8] global DOF numbers per element
# ndof   - total number of DOFs

function create_quad_mesh(nr, ntheta, ri = 5.0, ro = 10.0)
    nnodes = (nr + 1) * (ntheta + 1)    # number of nodes
    nelem = nr * ntheta                 # number of quadrilateral elements
    ndof = 2 * nnodes                   # 2 DOFs per node (ux, uy)

    # nodal coordinates with (polar -> cartesian)
    coords = zeros(nnodes, 2)
    dr = (ro - ri) / nr
    dtheta = (pi / 2) / ntheta

    for j in 0:ntheta
        theta = j * dtheta
        for i in 0:nr
            r = ri + i * dr
            n = j * (nr + 1) + i + 1

            coords[n, 1] = r * cos(theta)   # x1-coordinate
            coords[n, 2] = r * sin(theta)   # x2-coordinate
        end
    end

    # connectivity and EFT
    conn = zeros(Int, nelem, 4)     # node numbers [n1, n2, n3, n4]
    eft = zeros(Int, nelem, 8)      # global DOF numbers [ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]

    e = 0   # element counter
    for j in 0:ntheta-1
        for i in 0:nr-1
            # counter-clockwise node numbering
            n1 = j       * (nr + 1) + i       + 1       # lower left
            n2 = j       * (nr + 1) + (i + 1) + 1       # lower right
            n3 = (j + 1) * (nr + 1) + (i + 1) + 1       # upper right
            n4 = (j + 1) * (nr + 1) + i       + 1       # upper left

            e += 1
            conn[e, :] = [n1, n2, n3, n4]
            eft[e, :] = [2 * n1 - 1, 2 * n1,
                         2 * n2 - 1, 2 * n2,
                         2 * n3 - 1, 2 * n3,
                         2 * n4 - 1, 2 * n4]
        end 
    end

    return coords, conn, eft, ndof
end