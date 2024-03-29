from dolfin import *
import sympy as smp
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt

parameters["form_compiler"]["quadrature_degree"]=4


## Level set as a function to define the domain ==> not specifically needed
def phiFunc(x):
    return -(1.0/8.0 - (x[0]-0.5)**2 - (x[1]-0.5)**2)

## Solves the Poisson equation using phi-FEM without stabilisation terms
def solvePoissonStab(nCells=8, plotSol=False):

    ## The mesh can also be created like this
    big_mesh = UnitSquareMesh(nCells, nCells)
    domain = MeshFunction("size_t", big_mesh, big_mesh.topology().dim(), 0)
    domain.set_all(0)
    nbC = 0
    for c in cells(big_mesh):
        coords = c.get_vertex_coordinates()
        if phiFunc(coords[0:2]) < 0 or phiFunc(coords[2:4]) < 0 or phiFunc(coords[4:6]) < 0:
            # print("Domain:", c)
            nbC +=1
            domain[c] = 1
    mesh = SubMesh(big_mesh, domain, 1)

    ## Define function space
    V = FunctionSpace(mesh, 'CG', 1)

    ## Define level-set expression positive
    phi = Expression('-(1.0/8.0 - pow(x[0]-0.5, 2) - pow(x[1]-0.5, 2))',degree=1,domain=mesh)
    phi = interpolate(phi, V)

    ## Ghost penalty
    mesh.init(1,2)      ## enable connectivity ==> all facets connected to a given cell
    ghostCells = MeshFunction("size_t", mesh, mesh.topology().dim())
    ghostCells.set_all(0)
    ghostFacets = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    ghostFacets.set_all(0)

    nbGhost = 0
    for c in cells(mesh):           ## set ghost cells
        for f in facets(c):
            v1, v2 = vertices(f)
            if phi(v1.point().x(), v1.point().y())*phi(v2.point().x(), v2.point().y()) < 0:
                ghostCells[c] = 1
        if ghostCells[c] == 1:
            nbGhost +=1
            for f in facets(c):
                ghostFacets[f] = 1

    print("Total cells:", mesh.num_cells(), "   Ghost cells: ", nbGhost)

    dx = Measure("dx")(domain=mesh, subdomain_data=ghostCells)
    ds = Measure("ds")(domain=mesh)
    dS = Measure("dS")(domain=mesh, subdomain_data=ghostFacets)

 

    ## Define variational problem
    w = TrialFunction(V)        # equals u_D on the boundary
    v = TestFunction(V)         # equals 0 on the boundary

    # Computation of the Exact solution and exact source term
    x, y = smp.symbols('x1 x2')
    u_sympy = (0.125-(x-0.5)**2-(y-0.5)**2)*smp.exp(x)*smp.sin(2.0*smp.pi*y)       
    f_sympy = -smp.diff(smp.diff(u_sympy, x),x)-smp.diff(smp.diff(u_sympy, y),y)
    f = Expression(smp.ccode(f_sympy).replace('x1', 'x[0]').replace('x2', 'x[1]'),degree=1,domain=mesh)
    n = FacetNormal(mesh)
    sigma = 20.0
    h = CellDiameter(mesh)      
    h_avg =avg(h)       ## lorsqu'on est sur une facet ghost, prendre la moyenne des deux cellules concernées

    G = sigma*h_avg*dot(jump(grad(phi*w),n),jump(grad(phi*v),n))*dS(1) + sigma*(h**2)*div(grad(phi*w))*div(grad(phi*v))*dx(1)
    a = dot(grad(phi*w), grad(phi*v))*dx - dot(grad(phi*w), n)*phi*v*ds + G

    Grhs = sigma*(h**2)*f*div(grad(phi*v))*dx(1)
    l = f*phi*v*dx - Grhs

    ## Compute solution
    w_h = Function(V)
    solve(a == l, w_h)
    u_h = phi*w_h

    if plotSol:
        ## Plot mesh
        # p = plot(mesh)
        # theta = np.linspace(0,2*np.pi,100)
        # plt.plot(np.sqrt(1/8)*np.cos(theta)+0.5, np.sqrt(1/8)*np.sin(theta)+0.5)

        # # Save solution for Paraview
        # file = File("SolutionPVD.pvd")
        # file << u_h

        ## Export domain file
        # domain_file = File("domain.pvd")
        # domain_file << domain
            
        # Plot solution
        p = plot(u_h)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Poisson solution using Phi-FEM WITH stabilisation", fontsize="x-large", y=1.02)
        ## add a colorbar
        plt.colorbar(p)
        ## save image to disk
        plt.savefig("Solution.png")
        ## Show the plot, on a server

        plt.show()

    ## Compute max cell edge length
    hMax = mesh.hmax()

    ## Compute exact solution 
    u_exact = Expression(smp.ccode(u_sympy).replace('x1', 'x[0]').replace('x2', 'x[1]'),degree=4,domain=mesh)

    ## Compute errors in L2 and H1 norms
    errExactL2 = assemble(u_exact**2*dx(0))**0.5 
    error_L2 = assemble((u_exact-u_h)**2*dx(0))**0.5 / errExactL2
    
    errExactH1 = assemble((u_exact)**2*dx(0) + dot(grad(u_exact), grad(u_exact))*dx(0))**0.5
    error_H1 = assemble((u_exact-u_h)**2*dx(0) + dot(grad(u_exact-u_h), grad(u_exact-u_h))*dx(0))**0.5 / errExactH1

    return hMax, error_L2, error_H1



if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        solvePoissonStab(nCells=100, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## wil contain errors to plot
        dataToPlot = []

        for i in range(1,5):
            print('/////////////////////////////////////////////')
            N = int(10*2**i)
            # compute errors
            hMax, errL2, errH1 = solvePoissonStab(N, False)
            # Print errors
            print('N: ', N, ' hMax: {:.3f}'.format(hMax), 
                ' L2 error: {:.3f}'.format(errL2), ' H1 error: {:.3f}'.format(errH1))
            # add to data to be plotted
            dataToPlot.append([hMax, errL2, errH1])

        ## Extract columns to be plotted
        dataToPlot = np.array(dataToPlot)
        hMax = dataToPlot[:, 0]
        errL2 = dataToPlot[:, 1]
        errH1 = dataToPlot[:, 2]
        ## Linear regression to estimate the slopes
        regR2 = spst.linregress(np.log10(hMax), np.log10(errL2))
        regH1 = spst.linregress(np.log10(hMax), np.log10(errH1))

        ## Set plot style and plot graphs
        plt.style.use("seaborn")
        # plt.loglog(hMax,10*hMax)
        plt.loglog(hMax, errL2, "bo-", label = "$L^2$ slope = %.3f"%regR2.slope)
        plt.loglog(hMax, errH1, "ro-", label = "$H^1$ slope = %.3f"%regH1.slope)

        ## Set title and other labels
        plt.title("$L^2$ and $H^1$ error convergence in log mode", fontsize="xx-large", y=1.025)
        plt.xlabel("$h_{max}$")
        plt.ylabel("$\Vert u - u_h \Vert$")
        plt.legend(fontsize="large")

        ## Save the figure
        plt.savefig("EtudeDeCvg.png")

        ## Show the plot in web browser
        plt.show()
