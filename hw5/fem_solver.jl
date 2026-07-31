# Assembly, BCs, solver, VTK export

# CONTENTS:
# assemble_system           - global K and f from element contributions
# apply_traction_top!       - distributed load on upper edge (load case a)
# apply_point_load!         - point load at a node (load case b)
# apply_bcs!                - enforce displacement boundary conditions
# find_bottom_dofs          - DOFs on the bottom edge (uy = 0)
# find_lower_left_x_dof     - DOF at the lower left node (ux = 0)
# find_upper_right_node     – node at the upper-right corner
# compute_stresses          – element stresses
# rite_vtk_results          – ParaView output via WriteVTK.jl

using FEMLab

# ASSEMBLY
function assemble_system(coords, conn, eft, ndof, E, nu, h; b = [0.0, 0.0])
    nelem = size(conn, 1)

    # triplet lists for sparse assembly
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    f = zeros(ndof)

    has_body = norm(b) > 0.0
    
    for e in 1:nelem
        nodes = conn[e, :]
        xe = coords[nodes, 1]
        ye = coords[nodes, 2]
        dofs = eft[e, :]

        # element stiffness matrix [6x6]
        Ke = element_stiffness(xe, ye, E, nu, h)

        # collect triplets
        for i in 1:6, j in 1:6
            push!(I_idx, dofs[i])
            push!(J_idx, dofs[j])
            push!(V_vals, Ke[i, j])
        end

        # element body force vector (only if body force is present)
        if has_body
            fe = element_body_force(xe, ye, b, h)
            f[dofs] .+= fe
        end
    end

    K = sparse(I_idx, J_idx, V_vals, ndof, ndof)

    return K, f
end

# DISTRIBUTED TRACTION ON UPPER EDGE (load case a)
function apply_traction_top!(f, nx, ny, t_vec, h, L)
    tx, ty = t_vec
    dx = L / nx

    for i in 0:nx
        node = ny * (nx + 1) + i + 1    # node number on upper edge
        dof_x = 2 * node - 1
        dof_y = 2 * node

        # corner -> half segment weight, interior -> full segment weight
        weight = (i == 0 || i == nx) ? 0.5 : 1.0

        f[dof_x] += h * tx * dx * weight
        f[dof_y] += h * ty * dx * weight
    end
end

# POINT LOAD AT NODE (load case b)
# Load case b: point load P in x-direction at the upper right corner
function apply_point_load!(f, node, direction, value)
    dof = (direction == 1) ? 2 * node - 1 : 2 * node
    f[dof] += value
end

# ENFORCE DISPLACEMENT BOUNDARY CONDITIONS
function apply_bcs!(K, f, fixed_dofs)
    for d in fixed_dofs
        K[d, :] .= 0.0
        K[:, d] .= 0.0
        K[d, d]  = 1.0
        f[d]     = 0.0
    end
end

# HELPER FUNCTIONS: DOF and node lookup
# u_y = 0 for all nodes on the bottom edge (j=0)
function find_bottom_dofs(nx, ny)
    dofs = Int[]
    for i in 0:nx
        node = i + 1          # j=0 → n = 0*(nx+1)+i+1
        push!(dofs, 2 * node) # u_y DOF
    end
    return dofs
end

# u_x = 0 at the lower left node (i=0, j=0 → node 1)
function find_lower_left_x_dof()
    return 1   # DOF 1 = u_x of node 1
end

# upper right corner: i=nx, j=ny
function find_upper_right_node(nx, ny)
    return ny * (nx + 1) + nx + 1
end

# STRESS RECOVERY
function compute_stresses(coords, conn, eft, u, E, nu)
    nelem = size(conn, 1)
    C = plane_stress_C(E, nu)
    sigma = zeros(nelem, 3)

    for e in 1:nelem
        nodes = conn[e, :]
        xe = coords[nodes, 1]
        ye = coords[nodes, 2]
        dofs = eft[e, :]

        B, _ = B_matrix(xe, ye)
        eps = B * u[dofs]
        sigma[e, :] = C * eps
    end 

    return sigma
end 

# VTK OUTPUT FOR PARAVIEW
function write_vtk_results(filename, coords, conn, u, sigma)
    nnodes    = size(coords, 1)
    nelem = size(conn, 1)
 
    # Points as 3×N matrix (ParaView expects 3D coordinates)
    points = zeros(3, nnodes)
    points[1, :] = coords[:, 1]     # x
    points[2, :] = coords[:, 2]     # y
    # points[3, :] = 0.0  (already zero from initialisation)
 
    # Displacements: u is [u_x1, u_y1, u_x2, u_y2, ...] -> reorder
    ux = u[1:2:end]
    uy = u[2:2:end]
    disp = zeros(3, nnodes)
    disp[1, :] = ux
    disp[2, :] = uy
 
    # VTK_TRIANGLE = cell type 5 in the VTK specification
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, conn[e, :]) for e in 1:nelem]
 
    vtk_grid(filename, points, cells) do vtkfile
        vtkfile["displacement", VTKPointData()] = disp
        vtkfile["sigma_x", VTKCellData()] = sigma[:, 1]
        vtkfile["sigma_y", VTKCellData()] = sigma[:, 2]
        vtkfile["sigma_xy", VTKCellData()] = sigma[:, 3]
    end

     println("  → Written: $(filename).vtu")
 end

