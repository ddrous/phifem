# PROJET M2 CSMI 2020 - MIMESIS

<span style="font-size:1.3em;">**Simulation des tissus mous à l’aide de méthodes éléments finis non-conformes innovantes (Phi-FEM)**</span>

<br>

## **Objectif**
L'objectif du projet est d'appliquer la méthode Phi-FEM à l'équation d'élasticité linéaire. Nous commencerons par comparer les résultats Phi-Phem avec la méthode classique sur l'équation de Poisson. 

<br>

## **Delivrables**
- Le rapport: `docs/pdfs/report-v2.pdf`
- Le code:
    - Poisson avec la méthode FEM classique: `src/ClassiFEM/Poisson.py`
    - Elasticité linéaire avec la méthode FEM classique: `src/ClassiFEM/Elasticity2D.py`
    - Poisson avec la méthode Phi-FEM: `src/PhiFEM/Poisson.py`
    - Elasticité linéaire avec la méthode Phi-FEM: `src/PhiFEM/Elasticity2D.py`

**Remarque:** Modifier la variable `CvgStudy` dans chacun de ces scripts pour obtenir soit une simple simulation, soit une étude de convergence. 

<br>

## **Ressources**
- [FEniCS](https://fenicsproject.org/): pour la formulation variationnelle
- [Docker](https://www.docker.com/): pour l'installation de FEniCS
- [Sympy](https://www.sympy.org/en/index.html): pour le calcul des dérivés utilisées dans la formulation variationnelle