# RoughMob Method

## Overview
The RoughMob method is a new approach in molecular dynamics simulations, designed to predict the artificial acceleration of dynamics when one changes from an atomistic (AA) to a coarse-grained (CG) simulation model. As of now, this method has been proven sucessfull in accurately estimating the acceleration factors of self-diffusion and binary (mutual) diffusion for pure hydrocarbon liquids and their binary mixtures. 

## Methodology
RoughMob method's core concept lies in correlating the molecular roughness difference (a quantity derived from a numerical comparison of the AA and CG molecular surfaces) and roughness volumes (derived from the molecular roughness difference), with dynamic acceleration. This correlation allows for an accurate prediction of the acceleration factor upon coarse-graining, crucial in determining "true" diffusion coefficients in CG molecular dynamic simulations. The method has been tested and verified across multiple alkanes and aromatic molecules and alkane combinations. 

## Further Reading and References

For more detailed information on the RoughMob method and its applications, please refer to the following key references:

1. Meinel, M. K., & Müller-Plathe, F. (2020). Loss of Molecular Roughness upon Coarse-Graining Predicts the Artificially Accelerated Mobility of Coarse-Grained Molecular Simulation Models. *Journal of Chemical Theory and Computation*, 16(3), 1411-1419. [DOI: 10.1021/acs.jctc.9b00943](https://doi.org/10.1021/acs.jctc.9b00943)

2. Meinel, M. K., & Müller-Plathe, F. (2022). Roughness Volumes: An Improved RoughMob Concept for Predicting the Increase of Molecular Mobility upon Coarse-Graining. *Journal of Physical Chemistry B*, 126(20), 3737-3747. [DOI: 10.1021/acs.jpcb.2c00944](https://doi.org/10.1021/acs.jpcb.2c00944)

## Key Features
- **Predictive Accuracy**: Capable of predicting self-diffusion acceleration factors with minimal error, comparable to routine measurement accuracies.
- **Minimal Parameter Fitting**: Pure hydrocarbon systems and their mixtures utilize the same parameters.
- **Structural Preservation**: Employs Iterative Boltzmann Inversion (IBI) for maintaining structural integrity during coarse-graining, ensuring matching radial distribution functions (RDFs) and densities.
- **Minimal Atomistic Simulations**: Only short AA simulations for the generation of target RDFs of the pure hydrocarbons (and at equimolar mixture) are required - no additional AA simulations for geometrical calculations are needed. No AA simulations of mixtures at various concentrations are needed.

## Applicability
While so far primarily developed for unpolar molecules within a specific size range (6 to 13 carbon atoms) that can be coarse-grained into one CG bead per molecule, the RoughMob method shows promise for broader applications. Future enhancements aim to extend its applicability to e.g. polymers, and a wider chemical variety.

## Usage
The repository provides tools and examples for implementing the RoughMob method. Users can leverage this method to enhance the efficiency and accuracy of their simulations, particularly in the study of simple fluids and their mixtures.

## Contributions
We welcome any contributions and suggestions to improve the RoughMob method. Please refer to our contribution guidelines for more information or feel free to contact us. 

