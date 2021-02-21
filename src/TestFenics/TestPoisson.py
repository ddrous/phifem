from dolfin import *

# Using matplotlib to save solution data
# parameters ["plotting_backend"] = "matplotlib"
from matplotlib import pyplot


# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = File("poisson.pvd")
file << u

# save solution in XDMF format
# file = XDMFFile("poisson.xdmf")
# file.write(u , 0)


# Plot solution
p = plot(u)
# p = plot(u, interactive=True)

# set colormap
p.set_cmap("viridis")
p.set_clim(0.0, 1.0)
# add a title to the plot
pyplot.title("Poisson problem solution")
# add a colorbar
pyplot.colorbar(p) ;
# save image to disk
pyplot.savefig("poisson_sol.png")
