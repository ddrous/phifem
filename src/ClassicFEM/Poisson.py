from fenics import *
from dolfin import *
from mshr import *      ## For higher quality meshes
import sympy as smp
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt

parameters["form_compiler"]["quadrature_degree"]=4

## Solves the Poisson equation using phi-FEM without stabilisation terms
def solvePoissonStab(nCells=8, plotSol=False):

    ## The mesh can also be created like this
    domain = Circle(Point(0.5, 0.5), np.sqrt(1/8))
    mesh = generate_mesh(domain, nCells)

    # print(mesh.num_cells())
    # print(mesh.hmax())

    ## Define function space
    V = FunctionSpace(mesh, 'CG', 1)

    # phi = interpolate(phi, V)

    ## Comppute the second term f-------------------------------
    x0 = smp.symbols('x')
    x1 = smp.symbols('y')
    phiExprPy = -((1/8) - (x0-0.5)**2 - (x1-0.5)**2)
    uExExprPy = phiExprPy * (smp.exp(x0)*(smp.sin(2*smp.pi*x1)))

    lapUExPy = -smp.diff(uExExprPy, x0, 2) - smp.diff(uExExprPy, x1, 2)
    
    lapUExC = str(smp.ccode(lapUExPy))
    lapUExC = lapUExC.replace('exp', 'e')
    lapUExC = lapUExC.replace('x', 'x[0]')
    lapUExC = lapUExC.replace('y', 'x[1]')
    # lapUExC = lapUExC.replace('pi', 'M_PI')
    lapUExC = lapUExC.replace('e', 'exp')
    # print("Laplacian of U_exact: ", lapUExC)
    ##----------------------------------------------------------

    ##---------------------------------------------------------------
    uExExprC = str(smp.ccode(uExExprPy))
    # print("AS STRING", str(lapUExPy))
    uExExprC = uExExprC.replace('exp', 'e')
    uExExprC = uExExprC.replace('x', 'x[0]')
    uExExprC = uExExprC.replace('y', 'x[1]')
    # lapUExC = lapUExC.replace('pi', 'M_PI')
    uExExprC = uExExprC.replace('e', 'exp')
    # print("U_exact Expression: ", uExExprC)
    ##----------------------------------------------------------------

    ## Define boundary condition expression
    u_D = Expression(uExExprC, degree=2)

    ## Identify BC: we will reuse the same boundaries FEniCS alrady knows
    def sameBoundary(x, on_boundary):
        return on_boundary

    ## Create the bc to use when solving the variational problem
    bc = DirichletBC(V, u_D, sameBoundary)

    ## Define variational problem
    u = TrialFunction(V)        # equals u_D on the boundary
    v = TestFunction(V)         # equals 0 on the boundary

    f = Expression(lapUExC, element=V.ufl_element(), domain=mesh)
    # f = Expression('-(exp(x)*((...) - 2*sin(2*)) )', element=V.ufl_element(), domain=mesh)
    n = FacetNormal(mesh)
    sigma = 500
    h = CellDiameter(mesh)      

    a = dot(grad(u), grad(v))*dx - dot(grad(u), n)*v*ds

    l = f*v*dx

    ## Compute solution
    u_h = Function(V)
    solve(a == l, u_h, bc)

    u_exact = Expression(uExExprC, degree=4, domain=mesh)
    u_exact = interpolate(u_exact, V)

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
        p = plot(u_exact)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Poisson solution using Classic FEM", fontsize="x-large", y=1.02)
        ## add a colorbar
        plt.colorbar(p)
        ## save image to disk
        plt.savefig("Solution.png")

        # ## Show the plot, on a server
        plt.show()

    ## Compute max cell edge length
    hMax = mesh.hmax()

    ##---------------------------------------------------------------

    ## Compute errors in L2 and H1 norms
    errExactL2 = assemble(u_exact**2*dx)**0.5 
    error_L2 = assemble((u_exact-u_h)**2*dx)**0.5 / errExactL2

    errExactH1 = assemble((u_exact)**2*dx + dot(grad(u_exact), grad(u_exact))*dx)**0.5

    error_H1 = assemble((u_exact-u_h)**2*dx + dot(grad(u_exact-u_h), grad(u_exact-u_h))*dx)**0.5 / errExactH1
    
    return hMax, error_L2, error_H1


if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        solvePoissonStab(nCells=100, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## will contain errors to plot
        dataToPlot = []

        for i in range(1,5):
            # N = 2**i
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
