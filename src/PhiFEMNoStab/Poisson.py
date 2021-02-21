from fenics import *
from dolfin import *
from mshr import *      ## For higher quality meshes
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt

## Solves the Poisson equation using phi-FEM without stabilisation terms
def solvePoissonNoStab(nCells=8, plotSol=False):

    ## Create a unit circle with nCells along the diameter
    # domain = Circle(Point(0, 0), 1)
    # mesh = generate_mesh(domain, nCells)        ## No need to define a submesh here because the level-set fits it perfectly

    ## Level set as a function to define the domain ==> not specifically needed
    def phiFunc(x):
        return (1 - x[0]*x[0] - x[1]*x[1]) / 4

    ## The mesh can also be created like this
    big_domain = Rectangle(Point(-1, -1), Point(1, 1))
    big_mesh = generate_mesh(big_domain, int(nCells*sqrt(2)))
    
    domain = MeshFunction("size_t", big_mesh, big_mesh.topology().dim(), 0)
    domain.set_all(0)

    for c in cells(big_mesh):
        coords = c.get_vertex_coordinates()
        if phiFunc(coords[0:2]) > 0 or phiFunc(coords[2:4]) > 0 or phiFunc(coords[4:6]) > 0:
            domain[c] = 1

    mesh = SubMesh(big_mesh, domain, 1)


    ## Define function space
    V = FunctionSpace(mesh, 'P', 1)

    ## Define level-set expression positive in the UnitCircle domain
    phi = Expression('(1 - x[0]*x[0] - x[1]*x[1]) / 4', element=V.ufl_element(), domain=mesh)       


    ## Define the exact solution (formaly the boundary)
    # u_D = Expression('(1 - x[0]*x[0] - x[1]*x[1]) * sin(x[0]) * exp(x[1]) / 4', degree=2)
    # w_D = Expression('sin(x[0]) * exp(x[1])', degree=2)

    ## Identify BC: we will reuse the same boundaries FEniCS alrady knows
    # def sameBoundary(x, on_boundary):
        # return on_boundary

    ## Create the bc to use when solving the variational problem
    # bc = DirichletBC(V, u_D, sameBoundary)
    # bc = DirichletBC(V, w_D, sameBoundary)

    ## Define variational problem
    w = TrialFunction(V)        # equals u_D on the boundary
    v = TestFunction(V)         # equals 0 on the boundary
    f = Expression('exp(x[1]) * (x[0]*cos(x[0]) + (1+x[1])*sin(x[0]))', element=V.ufl_element(), domain=mesh)
    n = FacetNormal(mesh)
    a = dot(grad(phi*w), grad(phi*v))*dx - dot(grad(phi*w), n)*phi*v*ds
    l = f*phi*v*dx

    ## Compute solution
    w_h = Function(V)
    solve(a == l, w_h)
    # solve(a == l, w_h)
    u_h = phi*w_h

    if plotSol:
        ## Plot solution and mesh: doent work on Docker
        p = plot(u_h)
        # plot(mesh)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Poisson solution using Phi-FEM WITHOUT stabilisation", fontsize="x-large", y=1.02)
        ## add a colorbar
        plt.colorbar(p)
        ## save image to disk
        plt.savefig("Solution.png")
        ## Show the plot, on a server
        plt.show()

    # ## Compute max cell edge length
    # hMax = mesh.hmax()
    # ## Compute errors in L2 and H1 norms
    # u_exact = u_D
    # error_L2 = errornorm(u_exact, u_h, 'L2')
    # error_H1 = errornorm(u_exact, u_h, 'H1')

    # return hMax, error_L2, error_H1


if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        solvePoissonNoStab(nCells=128, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## wil contain errors to plot
        dataToPlot = []

        for i in range(2,7):
            N = 2**i
            # compute errors
            hMax, errL2, errH1 = solvePoissonNoStab(N, False)
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
        plt.loglog(hMax, errL2, "bo-", label = "$L^2$ slope = %.3f"%regR2.slope)
        plt.loglog(hMax, errH1, "ro-", label = "$H^1$ slope = %.3f"%regH1.slope)

        ## Set title and other labels
        plt.title("$L^2$ and $H^1$ error convergence in log mode", fontsize="xx-large", y=1.025)
        plt.xlabel("$h_{max}$")
        plt.ylabel("$\Vert u - u_h \Vert$")
        plt.legend(fontsize="large")

        ## Save the figure
        plt.savefig("EtudeDeCvg.png")

