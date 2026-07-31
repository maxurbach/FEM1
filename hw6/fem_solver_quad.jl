# DESCRIPTION:
# Global solver procedures for 2D plane stress analysis using Quad4 elements

# Handles global stiffness assembly into sparse matrix format, imposition of homogeneous
# and inhomogeneous boundary conditions, stress recovery at element enters, and export 
# of displacement and stress fields to VTK files (.vtu)

# CONTENTS:
# assemble_quad_system   - global matrix K and load vector f
# apply_curved_beam_bcs! - enforces boundary conditions for the curved beam problem
# compute_quad_stresses  - computes element stresses at element center (xi=0, eta=0)
# write_quad_vtk_results - writes displacements and stresses to VTK file for visualization in ParaView

using FEMLab

include("plane_stress_C.jl")
include("quad_element.jl")

# ASSEMBLY
function assemble_quad_system(coords, conn, eft, ndof, E, nu, t)
    nelem = size(conn, 1)
    
    # triplet lists for sparse assembly
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]
    f = zeros(ndof)

    for e in 1:nelem
        nodes = conn[e, :]
        xe    = coords[nodes, 1]
        ye    = coords[nodes, 2]
        dofs  = eft[e, :]

        # element stiffness matrix [8x8]
        Ke = quad4_element_stiffness(xe, ye, E, nu, t)

         # collect triplets
        for i in 1:8
            for j in 1:8
                push!(I_idx, dofs[i])
                push!(J_idx, dofs[j])
                push!(V_vals, Ke[i, j])
            end
        end
    end

    K = sparse(I_idx, J_idx, V_vals, ndof, ndof)
    return K, f
end

# BOUNDARY CONDITIONS FOR CURVED BEAM PROBLEM
function apply_curved_beam_bcs!(K, f, coords, nr, ntheta, u0 = 0.01)
    prescribed_dofs = Int[]
    prescribed_vals = Float64[]

    # Boundary Gamma_II (theta = 0 -> j = 0)
    for i in 0:nr
        node = i + 1 # nodes on theta = 0 axis
        dof_x2 = 2 * node
        
        push!(prescribed_dofs, dof_x2)
        push!(prescribed_vals, u0) # u2 = u0 = 0.01
    end

    # Boundary Gamma_I (theta = pi/2 -> j = ntheta)
    for i in 0:nr
        node = ntheta * (nr + 1) + i + 1 # nodes on theta = pi/2 axis
        dof_x1 = 2 * node - 1
        
        push!(prescribed_dofs, dof_x1)
        push!(prescribed_vals, 0.0) # u1 = 0
        
        # Pinned node at inner radius (r = ri, i = 0)
        if i == 0
            dof_x2 = 2 * node
            push!(prescribed_dofs, dof_x2)
            push!(prescribed_vals, 0.0) # u2 = 0
        end
    end

    # Modify RHS load vector f for inhomogeneous BCs (val != 0.0)
    for k in eachindex(prescribed_dofs)
        dof = prescribed_dofs[k]
        val = prescribed_vals[k]
        if val != 0.0
            f .-= K[:, dof] .* val
        end
    end

    # Apply Dirichlet boundary conditions (zero rows/cols, set diag = 1.0)
    for k in eachindex(prescribed_dofs)
        dof = prescribed_dofs[k]
        val = prescribed_vals[k]
        
        K[dof, :] .= 0.0
        K[:, dof] .= 0.0
        K[dof, dof] = 1.0
        f[dof] = val
    end
end

# STRESS RECOVERY
function compute_quad_stresses(coords, conn, eft, u, E, nu)
    nelem = size(conn, 1)
    C = plane_stress_C(E, nu)
    sigma = zeros(nelem, 3)

    for e in 1:nelem
        nodes = conn[e, :]
        xe = coords[nodes, 1]
        ye = coords[nodes, 2]
        dofs = eft[e, :]

        # Stress at element center (xi=0, eta=0)
        B, _ = quad4_B_and_detJ(xe, ye, 0.0, 0.0)
        eps = B * u[dofs]
        sigma[e, :] = C * eps
    end

    return sigma
end

# VTK OUTPUT FOR PARAVIEW
function write_quad_vtk_results(filename, coords, conn, u, sigma)
    nnodes = size(coords, 1)
    nelem = size(conn, 1)

    points = zeros(3, nnodes)
    points[1, :] = coords[:, 1]
    points[2, :] = coords[:, 2]

    ux = u[1:2:end]
    uy = u[2:2:end]
    disp = zeros(3, nnodes)
    disp[1, :] = ux
    disp[2, :] = uy

    # Cell type for bilinear Quad is VTKCellTypes.VTK_QUAD (cell type 9)
    cells = [MeshCell(VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL, conn[e, :]) for e in 1:nelem]

    vtk_grid(filename, points, cells) do vtkfile
        vtkfile["displacement", VTKPointData()] = disp
        vtkfile["sigma_x", VTKCellData()] = sigma[:, 1]
        vtkfile["sigma_y", VTKCellData()] = sigma[:, 2]
        vtkfile["sigma_xy", VTKCellData()] = sigma[:, 3]
    end

     println("  → Written: $(filename).vtu")
end
