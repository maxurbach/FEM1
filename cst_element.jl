# CST element routines

# The CST element has 3 nodes and 6 DOFs

# CONTENTS:
# plane_stress_C.    - elasticity matrix C for plane stress
# B_matrix           - strain-displacement matrix B (maps DOFs to strains)
# element_stiffness  - element stiffness matrix Ke = h * A * B' * C * B (plate thickness h)
# element_body_force - consistent nodal load vector for body forces

# ELASTICITY MATRIX (plane stress)
function plane_stress_C(E, nu)   
    f = E / (1 - nu^2)

    return f * [1.0  nu   0.0;
                nu   1.0  0.0;
                0.0  0.0  (1.0 - nu) / 2.0]
end

# STRAIN-DISPLACEMENT MATRIX (B_matrix)
# RETURNS: B [3x6], A (triangle area)
function B_matrix(xe, ye)
    x1, x2, x3 = xe
    y1, y2, y3 = ye

    # signed double area of the triangle (sign encodes node orientation)
    A2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    A = abs(A2) / 2.0

    if A < 1e-14
        error("Degenerate element with zero area detected")
    end

    b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2
    c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1

    B = (1.0 / A2) * [b1  0.0 b2  0.0 b3  0.0;
                      0.0 c1  0.0 c2  0.0 c3;
                      c1  b1  c2  b2  c3  b3]

    return B, A
end

# ELEMENT STIFFNESS MATRIX
# RETURNS: Ke [6x6]
function element_stiffness(xe, ye, E, nu, h)
    C = plane_stress_C(E, nu)
    B, A = B_matrix(xe, ye) 

    return h * A * B' * C * B
end

# CONSISTENT NODAL LOAD VECTOR FOR BODY FORCES
# The body force is distributed equally to all 3 nodes
# RETURNS: fe [6x1]
function element_body_force(xe, ye, b, h)
    bx, by = b
    _, A = B_matrix(xe, ye)

    return (h * A / 3.0) * [bx, by, bx, by, bx, by]
end