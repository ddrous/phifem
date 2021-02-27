# Phifem project

<span style="font-size:1.3em;">**Simulation of soft tissues using an innovative non-conforming Finite Element Methods (Phi-FEM)**</span>


## **Objectif**
The goal is to apply the novel Phi-fem technique to the elasticity equation, thus understanding how soft tissues (which are elastic), move. Before doing that, we will implement the Poisson equation and compare the results with the ones obtained using the classic FEM method.



## **Deliverables**
- The project report (in English): `docs/pdfs/report-v2.pdf`
- The codebase (in Python 3):
    - Poisson using classic FEM: `src/ClassiFEM/Poisson.py`
    - Linear elasticity using classic FEM: `src/ClassiFEM/Elasticity2D.py`
    - Poisson using Phi-FEM: `src/PhiFEM/Poisson.py`
    - Linear elasticity using Phi-FEM: `src/PhiFEM/Elasticity2D.py`

**Remarque:** Alter the variable`CvgStudy` in each of those scripts to perform whether a simple simulation, of a complete convergence study. 


## **Resources**
- [FEniCS](https://fenicsproject.org/): for solving the variational formulations
- [Docker](https://www.docker.com/): for installing FEniCS
- [Sympy](https://www.sympy.org/en/index.html): for computing the derivatives used in the variational formulations