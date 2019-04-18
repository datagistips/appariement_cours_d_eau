# Appariement des cours d'eau

Ce dossier comprend un script utilisé pour l'appariement de cours d'eau, en particulier ceux issus de deux lots distincts : des cours d'eau métiers issus de BDCarthage et ceux de la BDTOPO version 151.

Ces scripts ont été réalisés dans le cadre de l'identification des cours d'eau "Loi sur l'Eau" en 2015-2016, et dans le cadre d'un groupe de travail BDTOPAGE visant à la convergence de BDCarthage et de la BDTOPO. En effet, l'objet de l'unification de ces deux données était de bénéficier à la fois de la richesse sémantique de la première, et de la richesse géométrique de l'autre.

## Contenu
- `main.bat` : exécute le script `main.R`
- `main.R` : script "maître"
- `script.R` : exécuté par `main.R`, il contient les séries d'opérations visant à apparier les données
- `functions.R` : comprend le code des fonctions utilisées dans `script.R`

## Documentation
Le dossier `documentation` comprend le contenu de présentations réalisées notamment en séminaire, expliquant le procédé d'appariement basé sur la théorie des graphes

## Dossiers annexes
Le dossier `annexes_qualite` contient des scripts visant à qualifier a posteriori la qualité de l'appariement et cibler les endroits critiques où la distance entre le linéaire de référence et le linéaire apparié paraît excessive

Le dossier `annexes_tests` contient des scripts auxiliaires. L'opérationnalité de ces scripts n'est pas garantie, mais elle peut servir de source d'inspiration pour certains usages et calculs.

## Technologies
R a été utilisé pour ces calculs, en particulier la librairie igraph.

## Remarques
Le fichier `functions.R` comprend une fonction plutôt utile appelée `construireReseau` qui permet de créer un graphe depuis un réseau de linéaires et créant les noeuds et les arêtes nécessaires.
