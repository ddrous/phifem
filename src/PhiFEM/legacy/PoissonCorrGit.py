from __future__ import print_function
import numpy as np
from dolfin import *
import sympy
import matplotlib.pyplot as plt
#parameters['allow_extrapolation'] = True


# Number of iterations
init_Iter = 1
Iter = 1

# parameter of the ghost penalty
sigma = 20

# Polynome Pk
polV = 1
parameters["form_compiler"]["quadrature_degree"]=2*(polV+1)

# Ghost penalty
ghost = True

# plot the solution
Plot = True

# Compute the conditioning number
conditioning = False


def Omega(x,y):
	return -0.125+(x-0.5)**2+(y-0.5)**2<0


# Function used to write in the outputs files
def output_latex(f,A,B):
	for i in range(len(A)):
		f.write('(')
		f.write(str(A[i]))
		f.write(',')
		f.write(str(B[i]))
		f.write(')\n')
	f.write('\n')


# Computation of the Exact solution and exact source term
x, y = sympy.symbols('xx yy')
u1 = (0.125-(x-0.5)**2-(y-0.5)**2)*sympy.exp(x)*sympy.sin(2*sympy.pi*y)
f1 = -sympy.diff(sympy.diff(u1, x),x)-sympy.diff(sympy.diff(u1, y),y)


# Initialistion of the output
size_mesh_vec = np.zeros(Iter)
error_L2_vec = np.zeros(Iter)
error_H1_vec = np.zeros(Iter)
cond_vec = np.zeros(Iter)
for i in range(init_Iter-1,Iter):
	print('##################')
	print('## Iteration ',i+1,'##')
	print('##################')

	# Construction of the mesh
	N = int(10*2**((i)))
	mesh_macro = UnitSquareMesh(N,N)
	domains = MeshFunction("size_t", mesh_macro, mesh_macro.topology().dim())
	domains.set_all(0)
	for ind in range(mesh_macro.num_cells()):
		mycell = Cell(mesh_macro,ind)
		v1x,v1y,v2x,v2y,v3x,v3y = mycell.get_vertex_coordinates()
		if Omega(v1x,v1y) or Omega(v2x,v2y) or Omega(v3x,v3y):
			domains[ind] = 1
	mesh = SubMesh(mesh_macro, domains, 1)
	V = FunctionSpace(mesh,'CG',polV)

	# Construction of phi
	phi = Expression('-0.125+pow(x[0]-0.5,2)+pow(x[1]-0.5,2)',degree=polV,domain=mesh)
	phi = interpolate(phi,V)

	# Computation of the source term
	f_expr = Expression(sympy.ccode(f1).replace('xx', 'x[0]').replace('yy', 'x[1]'),degree=polV,domain=mesh)

	# Facets and cells where we apply the ghost penalty
	mesh.init(1,2)
	facet_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
	cell_ghost = MeshFunction("size_t", mesh, mesh.topology().dim())
	facet_ghost.set_all(0)
	cell_ghost.set_all(0)
	for mycell in cells(mesh):
		for myfacet in facets(mycell):
			v1, v2 = vertices(myfacet)
			if phi(v1.point().x(),v1.point().y())*phi(v2.point().x(),v2.point().y())<0:
				cell_ghost[mycell] = 1
				for myfacet2 in facets(mycell):
					facet_ghost[myfacet2] = 1

	# Initialize cell function for domains
	dx = Measure("dx")(domain = mesh,subdomain_data = cell_ghost)
	ds = Measure("ds")(domain = mesh)
	dS = Measure("dS")(domain = mesh,subdomain_data = facet_ghost)

	# Resolution
	n = FacetNormal(mesh)
	h = CellDiameter(mesh)
	u = TrialFunction(V)
	v = TestFunction(V)
	if ghost == False:
		a = inner(grad(phi*u),grad(phi*v))*dx - dot(inner(grad(phi*u),n),phi*v)*ds
		L = f_expr*v*phi*dx
	if ghost == True:
		a = inner(grad(phi*u),grad(phi*v))*dx - dot(inner(grad(phi*u),n),phi*v)*ds + sigma*avg(h)*dot(jump(grad(phi*u),n),jump(grad(phi*v),n))*dS(1)+sigma*h**2*inner(div(grad(phi*u)),div(grad(phi*v)))*dx(1)
		L = f_expr*v*phi*dx-sigma*h**2*inner(f_expr,div(grad(phi*v)))*dx(1)
	u_h = Function(V)
	solve(a == L, u_h,solver_parameters={'linear_solver': 'mumps'})
	sol = u_h*phi

	# Computation of the error
	u_expr = Expression(sympy.ccode(u1).replace('xx', 'x[0]').replace('yy', 'x[1]'),degree=4,domain=mesh)
	err_L2 = assemble((sol-u_expr)**2*dx(0))**0.5/assemble((u_expr)**2*dx(0))**0.5
	err_H1 = assemble((grad(sol-u_expr))**2*dx(0))**0.5/assemble((grad(u_expr))**2*dx(0))**0.5
	size_mesh_vec[i] = np.sqrt(2)/N
	error_L2_vec[i] = err_L2
	error_H1_vec[i] = err_H1
	print('h :',np.sqrt(2)/N)
	print('relative L2 error : ',err_L2)
	print('relative H1 error : ',err_H1)	
	if conditioning == True:
		A = np.matrix(assemble(a).array())
		ev, eV = np.linalg.eig(A)
		ev = abs(ev)
		cond = np.max(ev)/np.min(ev)
		cond_vec[i] = cond
		print("conditioning number x h^2",cond)
	print('')


# Print the output vectors
print('Vector h :',size_mesh_vec)
print('Vector relative L2 error : ',error_L2_vec)
print('Vector relative H1 error : ',error_H1_vec)
print("conditioning number",cond_vec)


#  Write the output file for latex
if ghost == False:
	f = open('output_no_ghost.txt','w')
if ghost == True:
	f = open('output_ghost2.txt','w')
f.write('relative L2 norm : \n')	
output_latex(f,size_mesh_vec,error_L2_vec)
f.write('relative H1 norm : \n')	
output_latex(f,size_mesh_vec,error_H1_vec)
if conditioning == True:
	f.write('conditioning number  : \n')	
	output_latex(f,size_mesh_vec,cond_vec)
f.close()


# Plot and save
if Plot == True:
	sol = project(u_h*phi,V)
	plot_sol = plot(sol)
	file = File('poisson.pvd')
	file << sol
	plt.savefig('myfig.png')