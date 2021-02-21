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

    L = 1.0      ## longeur de la poutre
    H = 1.0       ## hauteur

    rho = 1.0     ## f = -rho*9.81
    mu = 1.0      ## coeff de Lame
    lmda = 1.25 ## coeff de Lame lambda

    def phiFunc(x): ## Is this smooth enough ? YES !
        return -(1.0/8.0 - (x[0]-0.5)**2 - (x[1]-0.5)**2)

    ## Create a mesh macro for our elements
    aug = 0.00   ## This mesh is a bigger than the normal mesh
    big_mesh = RectangleMesh(Point(-aug,-aug), Point(L+aug,H+aug), int(L*nCells), nCells)
    
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

    ## Defnine Ghost cells and ghost facets
    mesh.init(1,2)      ## enable connectivities
    ghostCells = MeshFunction("size_t", mesh, mesh.topology().dim())
    ghostCells.set_all(0)
    ghostFacets = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    ghostFacets.set_all(0)

    nbGhost = 0
    for c in cells(mesh):           ## set ghost cells
        for f in facets(c):
            v1, v2 = vertices(f)
            v1Coord = [v1.point().x(), v1.point().y()]
            v2Coord = [v2.point().x(), v2.point().y()]
            if phiFunc(v1Coord) > 0 or phiFunc(v2Coord) > 0:
                ghostCells[c] = 1
    # for c in cells(mesh):           ## set ghost facets
        if ghostCells[c] == 1:
            nbGhost +=1
            # print("cellule: ", c)
            for f in facets(c):
                # print("Facet: ", f)
                ghostFacets[f] = 1
    print("Total cells:", nbC, "   Ghost cells: ", nbGhost)

    ## Defnie measures
    dx = Measure("dx")(domain=mesh, subdomain_data=ghostCells)
    ds = Measure("ds")(domain=mesh)
    dS = Measure("dS")(domain=mesh, subdomain_data=ghostFacets)

    ## Define function space
    V = VectorFunctionSpace(mesh, 'CG', 1, 2)

    phi = Expression('-(1.0/8.0 - pow(x[0]-0.5, 2) - pow(x[1]-0.5, 2))', degree=4, domain=mesh)       ## Phi is always scalar
    
    ## Define epsilon function for ease of use
    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    ## Define sigma function for ease of use
    def sigma(u):
        return lmda*div(u)*Identity(2) + 2*mu*epsilon(u)

    ## Define utility functions
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

    phiExpr = -((1.0/8.0) - (x0-0.5)**2 - (x1-0.5)**2)
    
    u1Expr = 2*x0 + smp.exp(x1)*smp.sin(x0)
    u2Expr = x0/2 + smp.cos(x0)-1

    # print("EXACT SOL: ", u1Expr, u2Expr)
    divU = smp.diff(u1Expr, x0, 1) + smp.diff(u2Expr, x1, 1)

    gradU1 = [smp.diff(u1Expr, x0, 1), smp.diff(u1Expr, x1, 1)]
    gradU2 = [smp.diff(u2Expr, x0, 1), smp.diff(u2Expr, x1, 1)]
    gradU = [gradU1, gradU2]

    Id = [[1.0,0.0],[0.0,1.0]]
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

    f = Expression((toCcode(f1), toCcode(f2)), element=V.ufl_element(), domain=mesh)

    u_exact = Expression((toCcode(u1Expr), toCcode(u2Expr)), element=V.ufl_element(), domain=mesh)

    g = u_exact + phi*Expression(('sin(x[0])', 'exp(x[1])'), degree=4, domain=mesh)

    ## Define variational problem
    w = TrialFunction(V)        # equals u_D on the boundary Gamma_Dphi
    v = TestFunction(V)         # equals 0 on the boundary

    uN = Constant((0, 0))

    n = FacetNormal(mesh)
    sigma_penalty = 20.0
    h = CellDiameter(mesh)      
    h_avg = (h('+') + h('-')) / 2. 


    G = sigma_penalty*h_avg*dot(jump(sigma(phi*w), n), jump(sigma(phi*v), n))*dS(1)
    G += sigma_penalty*(h**2)*dot(div(sigma(phi*w)), div(sigma(phi*v)))*dx(1)

    a = inner(sigma(phi*w), epsilon(phi*v))*dx - inner(dot(sigma(phi*w), n), phi*v)*ds + G

    Grhs = -sigma_penalty*(h**2)*(dot(f, div(sigma(phi*v))))*dx(1) - sigma_penalty*h_avg*dot(jump(sigma(g), n), jump(sigma(phi*v), n))*dS(1) - sigma_penalty*(h**2)*dot(div(sigma(g)), div(sigma(phi*v)))*dx(1)

    l = dot(f,phi*v)*dx - inner(sigma(g), epsilon(phi*v))*dx + inner(dot(sigma(g), n), phi*v)*ds + Grhs

    ## Compute solution
    w_h = Function(V)
    solve(a == l, w_h)

    u_h = phi * w_h + g
    u_h = project(u_h, V)
    u_exact = interpolate(u_exact, V)

    if plotSol:
        ## Plot mesh
        # p = plot(mesh)

        # # Save solution for Paraview
        file = File("SolutionPVD.pvd")
        file << project(u_h,V)

        ## Export domain file
        # domain_file = File("domain.pvd")
        # domain_file << domain

        # Plot solution
        p = plot(u_h)
        ## set colormap
        p.set_cmap("viridis")
        # p.set_clim(0.0, 1.0)  # to set vmax and vmin
        ## add a title to the plot
        plt.title("Elasticity solution using Phi FEM", fontsize="x-large", y=1.02)
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
        SolveElasticity(nCells=200, plotSol=True)

    ## We are solving multiple problems 
    else:
        ## wil contain errors to plot
        dataToPlot = []

        for i in range(2,6):
            # N = 2**i
            N = int(12*2**i)
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
