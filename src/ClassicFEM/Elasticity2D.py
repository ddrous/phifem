from fenics import *
from dolfin import *
from mshr import *      ## For higher quality meshes
import sympy as smp
import numpy as np
import scipy.stats as spst
from matplotlib import pyplot as plt

parameters["form_compiler"]["quadrature_degree"]=4

## Solves the Poisson equation using phi-FEM with stabilisation terms
def SolveElasticity(nCells=8, plotSol=False):

    L = 10.0      ## longeur de la poutre
    H = 1.0       ## hauteur
    rho = 1.0     ## f = -rho*9.81
    mu = 1.0      ## coeff de Lame
    lmda = 1.25 ## coeff de Lame lambda

    ## The mesh can also be created like this
    domain = Circle(Point(0.5, 0.5), np.sqrt(1.0/8.0))
    mesh = generate_mesh(domain, nCells)

    ## Define function space
    V = VectorFunctionSpace(mesh, 'P', 1, 2)

    ## Identify clamped endge
    def clampedBoundary(x, on_boundary):
        # return on_boundary and near(x[0],0.0)
        return on_boundary

    ## Defining clamped (Dirichlet) boundary
    notclamped = 0
    boundaries = MeshFunction("size_t", mesh, 1)
    boundaries.set_all(0);
    for f in facets(mesh):
        v1, v2 = vertices(f)
        v1Coord = [v1.point().x(), v1.point().y()]
        v2Coord = [v2.point().x(), v2.point().y()]
        if near(v1Coord[0],0.0) or near(v2Coord[0],0.0):
            notclamped+=1
            boundaries[f] = 1
    # print("Clamped facets:", notclamped, "  Total facets:", mesh.num_facets())
    ds = Measure("ds")(domain=mesh, subdomain_data=boundaries)

    ## Define epsilon function for ease of use
    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    ## Define sigma function for ease of use
    def sigma(u):
        return lmda*div(u)*Identity(2) + 2*mu*epsilon(u)


    def scalVec(scal, vec):
        return [vec[0]*scal, vec[1]*scal]

    def matVec(mat, vec):
        mat1 = mat[0][0]*vec[0], mat[0][1]*vec[1]
        mat2 = mat[1][0]*vec[0], mat[1][1]*vec[1]
        return [mat1, mat2]

    def matAdd(mat1, mat2):
        tmp1 = mat1[0][0]+mat2[0][0]
        tmp2 = mat1[0][1]+mat2[0][1]
        tmp3 = mat1[1][0]+mat2[1][0]
        tmp4 = mat1[1][1]+mat2[1][1]
        return [[tmp1, tmp2], [tmp3, tmp4]]

    def scalMat(scal, mat):
        tmp1 = scal*mat[0][0]
        tmp2 = scal*mat[0][1]
        tmp3 = scal*mat[1][0]
        tmp4 = scal*mat[1][1]
        return [[tmp1, tmp2], [tmp3, tmp4]]

    def transposeMat(mat):
        return [[mat[0][0], mat[1][0]], [mat[0][1], mat[1][1]]]

    ## calcul de u
    x0 = smp.symbols('x')
    x1 = smp.symbols('y')

    u1Expr = 2*x0 + smp.exp(x1)*smp.sin(x0)
    u2Expr = x0/2 + smp.cos(x0)-1

    divU = smp.diff(u1Expr, x0, 1) + smp.diff(u2Expr, x1, 1)

    gradU1 = [smp.diff(u1Expr, x0, 1), smp.diff(u1Expr, x1, 1)]
    gradU2 = [smp.diff(u2Expr, x0, 1), smp.diff(u2Expr, x1, 1)]
    gradU = [gradU1, gradU2]

    Id = [[1,0],[0,1]]
    sigma1 = scalMat(lmda*divU, Id)
    sigma2 = scalMat(mu, matAdd(gradU, transposeMat(gradU)))
    sigmaExpr = matAdd(sigma1, sigma2)

    f1 = -smp.diff(sigmaExpr[0][0], x0, 1) - smp.diff(sigmaExpr[0][1], x1, 1)
    f2 = -smp.diff(sigmaExpr[1][0], x0, 1) - smp.diff(sigmaExpr[1][1], x1, 1)

    def toCcode(pyExpr):
        cExpr = str(smp.ccode(pyExpr))
        cExpr = cExpr.replace('exp', 'e')
        cExpr = cExpr.replace('x', 'x[0]')
        cExpr = cExpr.replace('y', 'x[1]')
        cExpr = cExpr.replace('e', 'exp')
        return cExpr

    ## After computing the second term in PhiFEM, do not forget to come back here and compare
    f = Expression((toCcode(f1), toCcode(f2)), element=V.ufl_element(), domain=mesh)

    u_exact = Expression((toCcode(u1Expr), toCcode(u2Expr)), element=V.ufl_element(), domain=mesh)

    ## Create the bc to use when solving the variational problem
    ## Define boundary condition expression `u_D = g`
    u_D = u_exact       ## 3D vector
    bc = DirichletBC(V, u_D, clampedBoundary)

    ## Define variational problem
    u = TrialFunction(V)        # equals u_D on the boundary Gamma_D
    v = TestFunction(V)         # equals g on the boundary

    a = inner(sigma(u), epsilon(v))*dx
    l = dot(f,v)*dx

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
        p = plot(u_h)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Elasticity solution using Classic FEM", fontsize="x-large", y=1.02)
        ## add a colorbar
        plt.colorbar(p)
        ## save image to disk
        plt.savefig("Solution.png")

        # ## Show the plot, on a server
        plt.show()

    ## Compute max cell edge length
    hMax = mesh.hmax()
    # print(f1 ,"\n", f2)

    errExactL2 = assemble(inner(u_exact,u_exact)*dx)**0.5 
    error_L2 = assemble(inner((u_exact-u_h),(u_exact-u_h))*dx)**0.5 / errExactL2

    semiExactH1 = assemble(inner(grad(u_exact), grad(u_exact))*dx)**0.5
    semi_H1 = assemble(inner(grad(u_exact-u_h), grad(u_exact-u_h))*dx)**0.5 / semiExactH1

    ## Compute errors in L2 and H1 norms
    return hMax, error_L2, semi_H1


if __name__ == "__main__":
    
    ## What will you be doing?
    cvgStudy = False

    ## We are solving a single problem 
    if not cvgStudy:
        SolveElasticity(nCells=80, plotSol=True)

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
            print('N: ', N, ' hMax: {:.6f}'.format(hMax), 
                ' L2 error: {:.6f}'.format(errL2), ' H1 error: {:.6f}'.format(errH1))
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
        plt.loglog(hMax, errH1, "ro-", label = "Semi $H^1$ slope = %.3f"%regH1.slope)

        ## Set title and other labels
        plt.title("$L^2$ and $H^1$ error convergence in log mode", fontsize="xx-large", y=1.025)
        plt.xlabel("$h_{max}$")
        plt.ylabel("$\Vert u - u_h \Vert$")
        plt.legend(fontsize="large")

        ## Save the figure
        plt.savefig("EtudeDeCvg.png")

        ## Show the plot in web browser
        plt.show()
