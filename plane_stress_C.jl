# Elasticity matrix C for plane stress
# Copied from hw5/cst_element.jl

# ELASTICITY MATRIX (plane stress)
function plane_stress_C(E, nu)   
    f = E / (1 - nu^2)

    return f * [1.0  nu   0.0;
                nu   1.0  0.0;
                0.0  0.0  (1.0 - nu) / 2.0]
end