from fenics import *
from dolfin import *
from mshr import *      ## For higher quality meshes
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt


def solvePoisson(nCells=8, plotSol=False):

    ## Create a mesh: a unit square with nCells x nCells cells
    # mesh = UnitSquareMesh(nCells, nCells)

    ## Create a mesh: a unit circle with nCells along the diameter
    x_0 = 0.0
    y_0 = 0.0
    radius = 1.0
    domain = Circle(Point(x_0, y_0), radius)
    mesh = generate_mesh(domain, nCells)

    ## Define function space
    V = FunctionSpace(mesh, 'P', 1)

    ## Define boundary condition expression
    # u_D = Expression('(1 - x[0]*x[0] - x[1]*x[1])', degree=2)
    u_D = Expression('(1-x[0]*x[0]-x[1]*x[1]) * sin(x[0])*exp(x[1]) / 4', degree=2)

    ## Identify BC: we will reuse the same boundaries FEniCS alrady knows
    def sameBoundary(x, on_boundary):
        return on_boundary

    ## Identify BC: the solution will be u_D in the defines zones
    # tol = 1e-14
    # def zeroBoundary(x, on_boundary):
    #     return abs(x[0]) < tol or abs(x[1]) < tol or abs(x[0] - 1) < tol or abs(x[1] - 1) < tol

    ## Create the bc to use when solving the variational problem
    bc = DirichletBC(V, u_D, sameBoundary)

    ## Define variational problem
    u = TrialFunction(V)        # equals u_D on the boundary
    v = TestFunction(V)         # equals 0 on the boundary
    # f = Constant(1.0)
    f = Expression('exp(x[1]) * (x[0]*cos(x[0]) + (1+x[1])*sin(x[0]))', element=V.ufl_element(), domain=mesh)
    a = dot(grad(u), grad(v))*dx
    l = f*v*dx

    ## Compute solution
    u_h = Function(V)
    solve(a == l, u_h, bc)

    if plotSol:
        ## Plot solution and mesh: doent work on Docker
        p = plot(u_h)
        # plot(mesh)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Solution to the Poisson equation using classic FEM", fontsize="x-large", y=1.02)
        ## add a colorbar
        plt.colorbar(p)
        ## save image to disk
        plt.savefig("Solution.png")
        ## Show the plot, on a server
        plt.show()


    ## Compute max cell edge length
    hMax = mesh.hmax()
    ## Compute errors in L2 and H1 norms
    u_exact = u_D
    error_L2 = errornorm(u_exact, u_h, 'L2')
    error_H1 = errornorm(u_exact, u_h, 'H1')

    return hMax, error_L2, error_H1


if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        solvePoisson(nCells=128, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## wil contain errors to plot
        dataToPlot = []

        for i in range(2,7):
            N = 2**i
            # compute errors
            hMax, errL2, errH1 = solvePoisson(N, False)
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

