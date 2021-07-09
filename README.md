# Project Title : Assessment of dimensionality reduction methods for the detection of intermediate cancer phenotypes from 'omic' data

This repository is dedicaced to my master's degree internship project at the International Agency for Research on Cancer (IARC) with the Genetic Cancer Susceptibility Group.

Cancer classifications are a very important tool in cancer clinical treatments as it allows a precise tumor diagnostic and therefore the prescription of the most adequate and effective treatments available.
Recent breakthroughs in sequencing technology have lead to discovering that different cancer types of a same tissue do not exist in distinct classes but are rather linked to one an other in a spectrum of possible phenotypes. This discovery raises the concern of properly detecting rare cancers of intermediate phenotype as they have a higher risk of being miss classified or may require a very specific treatment.

In this context, this project goal is to assess how accurately cancer of intermediate phenotype can be detected using transcription data and dimensionality reduction methods (PCA and UMAP).

The project is devided into the following tasks:

1) Create a simulation script to simulate transcription data of :
- Two clusters of cancer samples that represent known classified cancers.
- Samples of rare cancers of intermediate phenotype (mixture of gene expression of the cluster samples)
[Simulation script (clusterSimulation.R)](scripts/clusterSimulation.R)

2) Create a script that produces a PCA & UMAP dimension reduction (DR), the DR projections of the simulated intermediate phenotype cancers of known mixture are compared to their expected projection, the difference between experimental and theoretical projections is measured using the R squared method.
This script allows to explore different UMAP parameters to assess how the average R square result is affect by those parameters.
The script outputs a boxplot of R square values for each parameter. 
[Boxplot script](scripts/simulationsBoxplot_PCA_UMAP.R)(simulationsBoxplot_PCA_UMAP.R)
Example of script output : 
[Simulation example](plot_examples/simulation_param_n_neighbors.png)(plot_examples/simulation_param_n_neighbors.png), 
[Boxplot example](plot_examples/boxplot_param_n_neighbors.png)(plot_examples/boxplot_param_n_neighbors.png)

3) Pushing further step 2, we can represent the average R square values of several UMAP parameter combinations by creating a grid of heatmaps. The goal in this step is to try and assess what parameters of UMAP allows the most accurate projection of the intermediate phenotype cancers.
[Heatmap script](scripts/simulationsHeatmap.R)(scripts/simulationsHeatmap.R), 
Example of heatmap grid : [Heatmap example](plot_examples/heatmap_grid.png)(plot_examples/heatmap_grid.png)

You will find further details about the project and references in the report (in french) : [Report](final_report.pdf)(final_report.pdf)


## Getting Started

This project is coded in R. You will need a computer with good memory capabilities and access to several R packages to run the analysis.
List of required packages are detailed in the source code.


## Author

* **Jérome Poizat** - *Initial work* - [jerome.poizat@laposte.net](mailto:jerome.poizat@laposte.net)

## Acknowledgments

The list of [contributors](presentation_support.pdf)(slide 36) who participated in this project.
