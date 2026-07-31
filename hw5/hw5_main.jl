# 2 element mesh with load case a, convergence study load case b

using FEMLab
using Printf

include("mesh.jl")
include("cst_element.jl")
include("fem_solver.jl")

# MATERIAL AND GEOMETRY PARAMETERS
E = 10000.0   # Young's modulus [N/m^2]
nu = 0.2      # Poission's ratio
h = 0.01      # plate thickness [m]
L = 1.0       # side length [m]
t_top = 100.0 # surface traction, load case a [N/m^2]
P = 10.0      # point load, load case b [N]
bvol = -10.0  # body force [N/m^3]

# 2-ELEMENT MESH (load case a)
nx_a = 1
ny_a = 1
coords_a, conn_a, eft_a, ndof_a = create_mesh(nx_a, ny_a, L)

K_a, f_a = assemble_system(coords_a, conn_a, eft_a, ndof_a, E, nu, h)
apply_traction_top!(f_a, nx_a, ny_a, [0.0, t_top], h, L)

fixed_a = unique(vcat(find_bottom_dofs(nx_a, ny_a), [find_lower_left_x_dof()]))
apply_bcs!(K_a, f_a, fixed_a)

u_a = K_a \ f_a
sigma_a = compute_stresses(coords_a, conn_a, eft_a, u_a, E, nu)

write_vtk_results("loadcase_a_2elem", coords_a, conn_a, u_a, sigma_a)

# CONVERGENCE STUDY (load case b)
for nelem in [32, 128, 512, 2048]
    nxy = Int(round(sqrt(nelem / 2)))
    nx_b = nxy
    ny_b = nxy
    coords_b, conn_b, eft_b, ndof_b, = create_mesh(nx_b, ny_b, L)

    K_b, f_b = assemble_system(coords_b, conn_b, eft_b, ndof_b, E, nu, h; b = [0.0, bvol])
    apply_point_load!(f_b, find_upper_right_node(nx_b, ny_b), 1, P)

    fixed_b = unique(vcat(find_bottom_dofs(nx_b, ny_b), [find_lower_left_x_dof()]))
    apply_bcs!(K_b, f_b, fixed_b)

    u_b = K_b \ f_b
    sigma_b = compute_stresses(coords_b, conn_b, eft_b, u_b, E, nu)

    actual = 2 * nx_b * ny_b

    write_vtk_results("loadcase_b_$(actual)elem", coords_b, conn_b, u_b, sigma_b)
end
