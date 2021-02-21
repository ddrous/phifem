
## **Journal de bord**

### <u>Semaine 1 (14/10/2020 - 21/10/2020)</u>

**Travail à faire** (par ordre de priorité):
1. Installer [FEniCS](https://fenicsproject.org/download/) en utilisant Docker (Installer la version la version la plus récente) (voir `Notes.md`)
2. Résoudre l'équation de Poisson en utilisant la méthode FEM classique dans FEniCS (conformement à [ce tutoriel](https://fenicsproject.org/pub/tutorial/html/._ftut1004.html), voir `src/poisson/Poisson.py` et `src/poisson/Solution.png`)
3. Effectuer l'étude de convergence en normes H<sup>1</sup> et L<sup>2</sup> (voir `src/poisson/Poisson.py` et `src/poisson/EtudeDeVcg.png`)
4. Lire les documents relatifs à la méthode Phi-FEM
    - Description de la méthode: https://hal.archives-ouvertes.fr/hal-02521111 
    - Etude du cas Neumann: https://arxiv.org/abs/2003.11733


**Questions:**
> Q: Comment **afficher** directement les plot depuis le conteneur sans sauvegarder avec matplotlib 
>
> R: forwarder le résultat à un port et utiliser son navigateur (voir DockerSeptup)?

> Q: Paraview demande les données de maillage qui ne sont pas visibles en dehors du conteneur (exporter le maillage aussi?) 
>
> R: Problème special de [Paraview 5.8.0](https://discourse.paraview.org/t/cannot-open-vtu-files-with-paraview-5-8/3759/7) --> upgrade to 5.8.1



### <u>Semaines 2 et 3 (21/10/2020 - 03/11/2020)</u>

**Travail à faire**:
1. En prenant la fonction level-set = parabole et solution_exacte = parabole*sin(x)*exp(y), appliquer la méthode Phi-Phem à l'equation de Poisson (sans les termes de stabilisation)
2.  Rajouter les termes de stabilisation 
    > Le terme de stabilisation G (Ghost-Penalty) represente la difference entre les gradient de deux cellules adjacentes proches du bord; on determine les aretes du bord à l'aide de la level-set lorsqu'elle change de signe en utilisant MeshFuntion dans FEniCS (e.g. MeshFunction("size_t", mesh, 1))
3.  Lire l'article



### <u>Semaine 4 (03/11/2020 - 10/11/2020)</u>

**Travail à faire**:
1.  Vérifier que les bonnes cellules sont sélectionnées (plotter un petit maillage sur cellui d'un cercle). Le cercle peut être obtenu comme le maillage d'un cercle tres fin.
2.  Faire le cas test du papier (penser à utiliser Sympy)
3.  Lire l'article


### <u>Semaine 5 (10/11/2020 - 17/11/2020)</u>

**Travail à faire**:
1. Détailler le rapport (CutFEM et XFEM, poutre de validation, références, etc.)
2. Corriger le program pour Poisson PhiFEM (division par des doubles, interpolation de u_exact)


### <u>Semaine 6 (17/11/2020 - 24/11/2020)</u>

**Travail à faire**:
1. Implémenter en FEniCS de l’équation d’élasticité avec FEM Classique
2. Réfléchir à la formulation variationnelle de PhiFEM pour l’équation d’élasticité
3. Ne pas oublier de rajouter l'estimation d'erreur à priori dans le rapport


### <u>Semaine 7 (24/11/2020 - 01/12/2020)</u>

**Travail à faire**:
1. Ecrire la formulation
2. Faire le cas test (Prendre g=0, dirichlet sur toutes les bordures pour appliquer PHEFEM, On modifie un tout petit peu g aux bord (p. 16) pour ne pas avoir la solution exacte direct !)
3. Ecrire le rapport (Un autre avantage de PhiFEM c'est que les intégrales sont toutes complètes. par de problèmes d'applatissement!)
4. Lire l'article (Mais pas urgent de l'include dans le rapport) !!!

*** Question ***
- Comment choisir la level-set ? `abs(x[0]/L + x[1]/H) + abs(x[0]/L - x[1]/H) - 2.0`?
- Dérivabilité de `phi`?

### <u>Semaine 8 (01/12/2020 - 08/12/2020)</u>

**Travail à faire**:
1. Finir l'élasticité en PhiFEM

### <u>Semaine 9 (08/12/2020 - 15/12/2020)</u>

**Travail à faire**:
1. Corriger le bug de convergence de la méthode PhiFEM
2. Rajouter les termes aux bords
