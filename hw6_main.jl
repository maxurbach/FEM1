# Main script for hw6 

using FEMLab
using Printf

include("mesh_quad.jl")
include("quad_element.jl")
include("fem_solver_quad.jl")

# MATERIAL AND GEOMETRY PARAMETERS
E = 10000.0   # Young's modulus
nu = 0.25     # Poisson's ratio
t = 0.1       # plate thickness
ri = 5.0      # inner radius
ro = 10.0     # outer radius
u0 = 0.01     # prescribed vertical displacement at Gamma_II

println("==================================================")
println(" HW6 (a): Verification on 2x2 Element Mesh ")
println("==================================================")

nr_verif = 2
ntheta_verif = 2
coords_v, conn_v, eft_v, ndof_v = create_quad_mesh(nr_verif, ntheta_verif, ri, ro)

K_v, f_v = assemble_quad_system(coords_v, conn_v, eft_v, ndof_v, E, nu, t)
apply_curved_beam_bcs!(K_v, f_v, coords_v, nr_verif, ntheta_verif, u0)

u_v = K_v \ f_v
sigma_v = compute_quad_stresses(coords_v, conn_v, eft_v, u_v, E, nu)

write_quad_vtk_results("curved_beam_2x2", coords_v, conn_v, u_v, sigma_v)


println("\n==================================================")
println(" HW6 (b): Mesh Convergence Study ")
println("==================================================")

# List of mesh configurations to analyze (nr, ntheta)
mesh_configs = [
    (4, 10),
    (8, 20),
    (16, 40),
    (32, 80)
]

for (nr, ntheta) in mesh_configs
    nelem = nr * ntheta
    println("Running mesh configuration: nr = $nr, ntheta = $ntheta ($nelem elements)...")

    # Build system and solve
    coords, conn, eft, ndof = create_quad_mesh(nr, ntheta, ri, ro)
    K, f = assemble_quad_system(coords, conn, eft, ndof, E, nu, t)
    apply_curved_beam_bcs!(K, f, coords, nr, ntheta, u0)

    u     = K \ f
    sigma = compute_quad_stresses(coords, conn, eft, u, E, nu)

    # Save VTK results
    filename = @sprintf("curved_beam_%dx%d", nr, ntheta)
    write_quad_vtk_results(filename, coords, conn, u, sigma)
end