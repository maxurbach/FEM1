using FEMLab
using Printf

## Conversion factors:
#  1 in   = 0.0254 m
#  1 lbf  = 4.44822 N
#  1 Mpsi = 6.895 GPa

## preprocessing (a)
# matrix of node coordinates
nodes = 
    [  0  0 ;   # node 1 [x, y] [in]
      30 20 ;   # node 2
      30  0 ;   # node 3
      60 30 ;   # node 4
      60  0 ;   # node 5
      90 30 ;   # node 6
      90  0 ;   # node 7
     120 20 ;   # node 8
     120  0 ;   # node 9
     150  0 ]   # node 10
nodes = nodes .* 0.0254     # element-wise multiplication to convert inches to meters

# matrix of element nodes
elements = 
    [ 1  2 ;    # element 1
      1  3 ;    # element 2
      2  3 ;    # element 3
      2  4 ;    # element 4
      3  4 ;    # element 5
      3  5 ;    # element 6
      4  5 ;    # element 7
      4  6 ;    # element 8
      5  6 ;    # element 9
      5  7 ;    # element 10
      6  7 ;    # element 11
      6  8 ;    # element 12
      6  9 ;    # element 13
      7  9 ;    # element 14
      8  9 ;    # element 15
      8 10 ;    # element 16
      9 10 ]    # element 17

# Material properties and cross-sectional areas (same for all elements)
E_val = 10e6 * 6895.0  # Young's modulus in Pa (converted from Mpsi)
A_val = 1 * 0.0254^2   # cross-sectional area in m^2 (converted from in^2)

# number of nodes
nnodes = size(nodes, 1)

# number of elements
nelem = size(elements, 1)

# number of degrees of freedom
ndof = 2 * nnodes 

# vectors of bar parameters
bar_youngs_modulus = fill(E_val, nelem)
bar_area = fill(A_val, nelem)

# element freedom table (4 DOFs per element)
# DOFs of element i = (2i-1, 2i) 
EFT = zeros(Int, 4, nelem)
for k = 1:nelem
    ni = elements[k, 1]
    nj = elements[k, 2]
    EFT[:, k] = [2*ni-1, 2*ni, 2*nj-1, 2*nj]
end


## assembly (b)
local_element_stiffness_matrix = (E, A, L) -> begin
    (E * A / L) *
    [ 1  0 -1  0 ;
      0  0  0  0 ;
     -1  0  1  0 ;
      0  0  0  0 ]
end

transformation_matrix = (phi) -> begin
    c = cos(phi); s = sin(phi)
    [ c  s  0  0 ;
     -s  c  0  0 ;
      0  0  c  s ;
      0  0 -s  c ]
end

# allocate arrays for sparse matrix assembly
I = zeros(Int, 16 * nelem)
J = zeros(Int, 16 * nelem)
V = zeros(Float64, 16 * nelem)

# contribution index
ind = 1;

# iterate through all elements
for k = 1:nelem

    E_k = bar_youngs_modulus[k]
    A_k = bar_area[k]

    # node positions
    ni = elements[k, 1]
    nj = elements[k, 2]

    # compute metrics
    xi = nodes[ni, :]   # [x_i, y_i]
    xj = nodes[nj, :]   # [x_j, y_j]

    t = xj - xi                        # element vector
    L = norm(t)                        # element length
    phi = angle(complex(t[1], t[2]))   # orientation angle


    # compute local element stiffness matrix
    ke_local = local_element_stiffness_matrix(E_k, A_k, L)
    
    # compute transformation/globalization matrix
    T_mat = transformation_matrix(phi)

    # compute globalized element stiffness matrix
    Ke = T_mat' * ke_local * T_mat      # globalized 4x4 element stiffness matrix

    # collect element contributions
    for j = 1:4
        for i = 1:4
            V[ind] = Ke[i, j];
            I[ind] = EFT[i, k];
            J[ind] = EFT[j, k];

            # increment contribution index
            ind = ind + 1;
        end
    end
end

# define sparse stiffness matrix from (row, col, val) triplets
K = sparse(I, J, V, ndof, ndof)


## solution

# define load vector (c)
F = zeros(Float64, ndof)

# Node 3 (DOFs 5 and 6): Fy = -500 lbf
F[6] = -500

# Node 4 (DOFs 7 and 8): Fy = -1000 lbf
F[8] = -1000

# Node 6 (DOFs 11 and 12): Fx = 500 lbf, Fy = -1000 lbf
F[11] = 500
F[12] = -1000

# Node 9 (DOFs 17 and 18): Fy = -500 lbf
F[18] = -500

F = F .* 4.44822     # element-wise multiplication to convert lbf to Newtons


# apply homogeneous boundary conditions by modification (d)
# Node 1 (pin): Ux = 0, Uy = 0 (DOFs 1 and 2)
# Node 10 (roller): Uy = 0 (DOF 20)
constrained_dofs = [1, 2, 20]

K = Matrix(K)
for dof in constrained_dofs
    K[dof, :] .= 0.0
    K[:, dof] .= 0.0
    K[dof, dof] = 1.0
    F[dof] = 0.0
end
K = sparse(K)

# solve for displacements (e)
U = K \ F

println("=" ^ 55)
println("Nodal Displacements")
println("=" ^ 55)
println(@sprintf("%-6s  %12s  %12s  %12s", "Node", "Ux [m]", "Uy [m]", "|U| [m]"))
println("-" ^ 55)
 
max_disp = 0.0
max_node = 0
for i = 1:nnodes
    ux = U[2*i-1]
    uy = U[2*i]
    u_total = sqrt(ux^2 + uy^2)
    println(@sprintf("%-6d  %12.6e  %12.6e  %12.6e", i, ux, uy, u_total))
    if u_total > max_disp
        max_disp = u_total
        max_node = i
    end
end
println("=" ^ 55)
println(@sprintf("Maximum total displacement: %.6e m at node %d", max_disp, max_node))