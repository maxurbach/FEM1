using FEMLab

# Parameters
L = 6.0                # Length of the rod [m]
A = 10.0               # Cross-sectional area [m^2]
k = 5.0                # Thermal conductivity [J/(ms°C)]
T_hat = 5.0            # Dirichlet boundary condition at x=0 [°C]
q_hat = 15.0           # Neumann boundary condition at x=L [J/sm^2]
Q = 100.0              # Heat supply over the domain [J/sm]
nelem = 4              # Number of elements
nnodes = nelem + 1     # Number of nodes
l = L / nelem          # Length of each element (1.5 m)

# Local stiffness matrix and local load vector
ke = (k * A / l) * [1 -1; -1 1]
fe = (Q * l / 2) * [1.0; 1.0]

# Assembly
K = zeros(nnodes, nnodes)
f = zeros(nnodes)

# EFT
for k = 1:nelem
    dofs = [k, k + 1]
    K[dofs, dofs] += ke
    f[dofs]       += fe
end

# Neumann boundary condition at x = L
f[end] += A * q_hat

# Dirichlet boundary condition at x = 0 (non-homogeneous)
# Modification of the load vector and the stiffness matrix
f[2:end] -= K[2:end, 1] .* T_hat
# Replace the first row/column
K[1, :] .= 0.0
K[:, 1] .= 0.0
K[1, 1]  = 1.0
f[1]     = T_hat

# Solution
T = K \ f

# Exact solution (analytical)
x_fine = range(0, L, length=300)
T_exact = T_hat .+ (q_hat / (k) .+ Q*L/(k*A)) .* x_fine .- Q/(2*k*A) .* x_fine.^2
dT_exact = (q_hat/(k) + Q*L/(k*A)) .- Q/(k*A) .* x_fine

# FE solution (piecewise linear interpolation)
x_nodes = range(0, L, length=nnodes)
T_fe = T  # Node values

# Temperature gradient (piecewise constant per element)
x_elem_mid = [(e - 0.5) * l for e in 1:nelem]
dT_fe = [(T[e+1] - T[e]) / l for e in 1:nelem]

# Plots
p1 = plot(x_fine, T_exact, label="T exact", lw=2, color=:blue, xlabel="x [m]",
          ylabel="T [°C]", title="Temperature Distribution")
plot!(p1, x_nodes, T_fe, label="T FE", lw=2, color=:red,
      marker=:circle, linestyle=:dash)

p2 = plot(x_fine, dT_exact, label="dT/dx exact", lw=2, color=:blue, xlabel="x [m]",
          ylabel="dT/dx [°C/m]", title="Temperature Gradient")
scatter!(p2, x_elem_mid, dT_fe, label="dT/dx FE",
         color=:red, marker=:square, markersize=7)

plot(p1, p2, layout=(2,1), size=(700, 600))