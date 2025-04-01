# Analyses Multiblocks Interactives

Analyses Multiblocks Interactives est une application interactive développée en R avec Shiny, dédiée à l'analyse intégrée de données omiques. Elle permet aux chercheurs de travailler simultanément sur des jeux de données protéomiques, métabolomiques et de metadata pour identifier des tendances et des relations via des techniques d'analyse multivariée.

## Table des matières

- [Fonctionnalités](#fonctionnalités)
- [Installation](#installation)
- [Utilisation](#utilisation)
- [Contribution](#contribution)
- [Licence](#licence)
- [Auteur](#auteur)

## Fonctionnalités

- **Import de données :** Chargez facilement vos fichiers CSV pour les données protéomiques, métabolomiques et les metadata.
- **Prétraitement :** Nettoyage et normalisation des données pour une analyse optimale.
- **Analyses multivariées :**
  - **PCA** pour l'exploration des structures de données.
  - **PLS‑DA** et **OPLS** pour la classification et la discrimination.
  - **Block PLS‑DA** pour l'analyse conjointe de plusieurs jeux de données.
  - **PLS Regression** et **PLS Canonical** pour l'analyse de régression.
  - **ComDim** et **Consensus OPLS** pour des approches complémentaires.
- **Visualisations interactives :** Graphiques dynamiques réalisés avec ggplot2 et plotly.
- **Exportation :** Export des résultats d'analyse en fichiers Excel pour une utilisation ultérieure.

## Installation

1. **Cloner le dépôt :**

   ```bash
   git clone [<URL_DU_REPO>](https://github.com/axellebourez/AnalysesMultiblockInteractives)
   Ouvrir le projet dans RStudio.

2. **Installer les packages :**
Dans R, exécutez le code suivant :

r
Copier
install.packages(c("shiny", "shinydashboard", "data.table", "ggplot2",
                   "reshape2", "ropls", "mixOmics", "MBAnalysis", "dplyr",
                   "plotly", "ggrepel", "openxlsx", "ConsensusOPLS",
                   "RVAideMemoire", "abind", "visNetwork", "igraph"))
                   
2. **Lancer l'application :**
Dans RStudio, ouvrez le script principal et exécutez :

r
Copier
shiny::runApp("App.R")


## Utilisation 

Chargement des données :
Utilisez les champs dédiés dans l'interface pour importer vos fichiers CSV.

Sélection de l'analyse :
Parcourez les différents onglets pour lancer l'analyse souhaitée (PCA, PLS‑DA, OPLS, etc.).

Interactivité :
Explorez les graphiques interactifs et ajustez les paramètres selon vos besoins.

Exportation des résultats :
Téléchargez vos résultats au format Excel pour les partager ou les conserver.

## Contributions
Les contributions et suggestions sont les bienvenues !
Pour signaler un problème ou proposer une amélioration, veuillez ouvrir une issue ou soumettre une pull request sur GitHub.

## Licence
Ce projet est distribué sous licence MIT.

## Auteur
Axelle Bourez
