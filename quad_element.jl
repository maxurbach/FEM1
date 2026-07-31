# DESCRIPTION:
# Provides element routines for 4-node isoparametric quadrilateral elements (Quad4)

# Includes bilinear shape functions, isoparametric mapping, Jacobian evaluation, 
# strain-displacement matrix B, and numerical integration via 2x2 Gauss-Legendre quadrature

# CONTENTS:
# quad4_shape_functions    - the 4 shape functions and derivatives in natural coordinates
# quad4_shape_derivatives  - derivatives of the shape functions
# quad4_B_and_detJ         - strain-displacement matrix B and Jacobian determinant
# quad4_element_stiffness  - computes the 8x8 element stiffness matrix Ke using 2x2 Gauss quadrature

using FEMLab

# SHAPE FUNCTIONS N_i(xi, eta)
# RETURNS: N [4x1]
function quad4_shape_functions(xi, eta)
    N = 0.25 * [(1.0 - xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 + eta),
                (1.0 - xi) * (1.0 + eta)]
    return N
end

# DERIVATIVES dN_i / dxi and dN_i / deta
# RETURNS: dN [4x2] (columns: dN/dxi, dN/deta)
function quad4_shape_derivatives(xi, eta)
    dN = 0.25 * [-(1.0 - eta)  -(1.0 - xi);
                  (1.0 - eta)  -(1.0 + xi);
                  (1.0 + eta)   (1.0 + xi);
                 -(1.0 + eta)   (1.0 - xi)]
    return dN
end

# STRAIN-DISPLACEMENT MATRIX B and det(J)
# RETURNS: B [3x8], detJ
function quad4_B_and_detJ(xe, ye, xi, eta)
    dN = quad4_shape_derivatives(xi, eta)

    # Jacobian matrix J = [dx/dxi  dy/dxi; dx/deta  dy/deta] [2x2]
    Xelem = [xe ye]
    J = Xelem' * dN
    detJ = det(J)

    if detJ < 1e-14
        error("Degenerate element with zero Jacobian detected")
    end

    # Derivatives with respect to physical coordinates (dN/dx, dN/dy)
    invJ = inv(J)
    dN_phys = dN * invJ'

    # Assemble B-matrix [3x8]
    B = zeros(3, 8)
    for i in 1:4
        dNdx = dN_phys[i, 1]
        dNdy = dN_phys[i, 2]

        B[1, 2 * i -1]  = dNdx
        B[2, 2 * i]     = dNdy
        B[3, 2 * i - 1] = dNdy
        B[3, 2 * i]     = dNdx
    end

    return B, detJ
end

# ELEMENT STIFFNESS MATRIX
# RETURNS: Ke [8x8]
function quad4_element_stiffness(xe, ye, E, nu, t)
    C = plane_stress_C(E, nu)
    Ke = zeros(8, 8)

    # 2x2 Gauss integration points (+/- 1/sqrt(3))
    gp = 1.0 / sqrt(3.0)
    gauss_points = [-gp, gp]

    for xi in gauss_points
        for eta in gauss_points
            B, detJ = quad4_B_and_detJ(xe, ye, xi, eta)
            Ke .+= t * (B' * C * B) * detJ
        end
    end

    return Ke
end