## Étapes de mise en place de Docker pour utiliser FEniCS

Afin de faire fonctionner Fenics sur un des exemples de ce repository:
1. Installer Docker
2. Cloner ce repository 
3. Télécharger l'image Docker: `docker pull quay.io/fenicsproject/stable:latest`
4. Se placer à la racine du repo nouvellement crée
5. Lancer le conteneur: `docker run -ti -p 127.0.0.1:8000:8000 -v ${pwd}:"/home/fenics/shared" -w /home/fenics/shared --name fenics quay.io/fenicsproject/stable` (`$(pwd)` sous Linux). Notre répertoire (racine de ce repo) est alors monté dans le dossier `/home/fenics/shared` du conteneur.
6. Utiliser ces commandes par la suite si nécessaire (voir [ce tuto](https://fenics-containers.readthedocs.io/en/latest/introduction.html)): 
>```
>docker stop fenics-container   
>docker start fenics-container   
>docker exec -ti -u fenics fenics-container /bin/bash -l
>```
**Remarque:** Sous Visual Studio Code, une image FEniCS contenant en plus `scipy` et `matplotlib` a été configurée, et donc juste les étapes 1. et 2. précédentes sont nécessaires pour faire fonctionner les exemples (voir dossier `.devcontenair`).
**Attention:** Lors de l'utilisation de la fonction plot, Fenics demandera d'utiliser l'adresse `http://0.0.0.0:8000` sur le navigateur host. C'est inexact, il faudra utiliser `127.0.0.1:8000` comme indiqué dans la commande `docker run`.


## Utilisation de MeshFunction pour renumeroter les elements FEniCS

- >facet- > arete   
    >dx sur les cellules,  
    >ds sur les facet a l'exterieur (la bordure du domaine)   
    >dS sur les facet interieures au domaine

- identtifier les facet -> indentifier les vertex -> tester level set

    >MeshFunction retourne un fichier comme le maillage:
    >Pour les termes de stabilisation, il faut faire:
    >- 1 meshfuntion pour les facet
    >- 1 meshfuntion pour les cellules
    >
    >et les parcourir par une boucle imbriquée.

- To create mesh functions over the cell facets 
    ```
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    ```

- 
    ```
    for cell in cells(mesh): 
        c00[cell] = 1

    v1,v2=vertices(facet)
    phi(v1.x(),v1.())*phi(v2.x(),v2.y())
    ```

- Pour la stabilisation, des que phi est negatif, on garde tout ses points de la celulle, on obbtient un sous maillage
    ```
    submesh = SubMesh(mesh, new_region, 1) pour creer un maillage 
    ```
LIEN:
POUR LES SOUS-DOMAINES: https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/subdomains-poisson/python/documentation.html

POUR LA NOTATION UFL: https://fenics.readthedocs.io/projects/ufl/en/latest/manual/form_language.html

JUMP AND AVERAGE FENICS: https://fenicsproject.org/pub/course/lectures/2015-08-logg-beijing/lecture_10_discontinuous_galerkin.pdf
