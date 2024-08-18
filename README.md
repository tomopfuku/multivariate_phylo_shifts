# multivariate_phylo_shifts


## Description

This is a simple method for detecting shifts in the structure of continuous trait correlations along phylogenies. Details about the model are described by Parins-Fukuchi et al. (preprint available here: https://ecoevorxiv.org/repository/view/4645/). The general goal is to identify nodes on a phylogeny where the correlation matrix between continuous characters shifts, with the idea that such shifts might explain changes in the functional, genetic, and/or environmental dependencies between traits. 

Traits are modeled under Brownian motion and shift points in the covariance matrix are identified using an adaptation of the algorithm described by Smith et al. 2023.



## Tutorial

The method is straightforward to run and requires only a tree and trait data. The tree should be represented in the newick format and the traits can be saved as a comma-separated values file formatted similarly to the example file here. Note that this method is still somewhat experimental and as such, I do not recommend attempting to use it without getting in touch first. 

To run the algorithm, use the `mod_evo.py` script located in the `code` folder. It can be run like so:

```
cd vitaceae
python ../code/mod_evo.py vitaceae.dated.macroprune.tre vitaceae.Leaves.impute.csv 10

```

The first two arguments are self-explanatory: the tree and the trait file. The last one is the minimum subclade size to consider. In all likelihood, the model will be underpowered to reasonablly reconstruct correlation matrices for clades that are too small and so we tell the search algorithm to ignore any shifts at a node leading to a clade smaller than our minimum threshold. For this dataset, we set the minimum clade size to 10 species. Note that two example files are given. One contains missing data, and the other has had missing measurements imputed using a phylogenetic imputation model, similar to Goolsby et al. (2017). Running each yields results that only slightly differ. This is reassuring. You may experiment with both, but in my experience, if the results differ significantly, that alone should be reason to take pause. 


**References**

Goolsby, E.W., Bruggeman, J. and Ané, C., 2017. Rphylopars: fast multivariate phylogenetic comparative methods for missing data and within‐species variation. Methods in Ecology and Evolution, 8(1), pp.22-27.

Parins-Fukuchi, C.T., Stull, G.W., Wen, J. and Beaulieu, J., 2022. Temperate-tropical transitions are linked with shifts in the structure of evolutionary integration in Vitaceae leaves.

Smith, S.A., Walker‐Hale, N. and Parins‐Fukuchi, C.T., 2023. Compositional shifts associated with major evolutionary transitions in plants. New Phytologist, 239(6), pp.2404-2415.

 
