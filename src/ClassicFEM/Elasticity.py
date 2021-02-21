from fenics import *
from dolfin import *
from mshr import *      ## For higher quality meshes
import sympy as smp
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt

parameters["form_compiler"]["quadrature_degree"]=4

## Solves the Poisson equation using phi-FEM without stabilisation terms
def SolveElasticity(nCells=8, plotSol=False):

    L = 10      ## longeur de la poutre
    W = 1       ## section carré
    rho = 1     ## f = -rho*9.81
    mu = 1      ## coeff de Lame
    lmda = 1.25 ## coeff de Lame lambda

    ## The mesh can also be created like this
    mesh = BoxMesh(Point(0,0,0), Point(L,W,W), 10*nCells, nCells, nCells)

    ## Define function space
    V = VectorFunctionSpace(mesh, 'P', 1)

    ## Define boundary condition expression
    u_D = Constant((0,0,0))       ## 3D vector

    ## Identify clamped endge
    def clampedBoundary(x, on_boundary):
        return on_boundary and near(x[0],0.0)

    ## Create the bc to use when solving the variational problem
    bc = DirichletBC(V, u_D, clampedBoundary)

    # ## Class for defning Neumann condition
    # class NotClamped(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[0] > tol

    # ## Define Neumann boundary conditions
    # neumann = NotClamped()
    # boundaries = MeshFunction("size_t", mesh, 2)
    # boundaries.set_all(0)
    # neumann.mark(boundaries, 1)
    # ds = Measure("ds")(domain=mesh, subdomain_data=neumann)

    ## Defining Neumann boundary
    notclamped = 0
    boundaries = MeshFunction("size_t", mesh, 2)
    boundaries.set_all(0);
    for f in facets(mesh):
        v1, v2, v3 = vertices(f)
        v1Coord = [v1.point().x(), v1.point().y(), v1.point().z()]
        v2Coord = [v2.point().x(), v2.point().y(), v2.point().z()]
        v3Coord = [v3.point().x(), v3.point().y(), v3.point().z()]
        if near(v1Coord[0],0.0) or near(v2Coord[0],0.0) or near(v3Coord[0],0.0):
            notclamped+=1
            boundaries[f] = 1
    print("Clamped facets:", notclamped, "  Total facets:", mesh.num_facets())
    ds = Measure("ds")(domain=mesh, subdomain_data=boundaries)

    ## Define epsilon function for ease of use
    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    ## Define sigma function for ease of use
    def sigma(u):
        return lmda*div(u)*Identity(3) + 2*mu*epsilon(u)

    ## Define variational problem
    u = TrialFunction(V)        # equals u_D on the boundary Gamma_D
    v = TestFunction(V)         # equals 0 on the boundary

    f = Constant((0,0,-rho*9.81))
    g = Constant((0,0,-1))

    a = inner(sigma(u), epsilon(v))*dx
    l = dot(f,v)*dx + dot(g,v)*ds(1)

    ## Compute solution
    u_h = Function(V)
    solve(a == l, u_h, bc)

    if plotSol:
        ## Plot mesh
        # p = plot(mesh)

        # # Save solution for Paraview
        file = File("SolutionPVD.pvd")
        file << u_h

        ## Export domain file
        # domain_file = File("domain.pvd")
        # domain_file << domain
            
        # Plot solution
        # p = plot(u_h)
        # ## set colormap
        # p.set_cmap("viridis")
        # # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        # ## add a title to the plot
        # plt.title("Poisson solution using Classic FEM", fontsize="x-large", y=1.02)
        # ## add a colorbar
        # plt.colorbar(p)
        # ## save image to disk
        # plt.savefig("Solution.png")

        # ## Show the plot, on a server
        # plt.show()

    ## Compute max cell edge length
    hMax = mesh.hmax()

    ## Compute errors in L2 and H1 norms
    return hMax, None, None


if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        SolveElasticity(nCells=10, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## wil contain errors to plot
        dataToPlot = []

        for i in range(1,5):
            # N = 2**i
            N = int(10*2**i)
            # compute errors
            hMax, errL2, errH1 = SolveElasticity(N, False)
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

        ## Show the plot in web browser
        plt.show()
