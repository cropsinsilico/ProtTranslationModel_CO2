# Protein translation model

## Model Description

This model essentially predicts the steady state protein concentration of selected enzymes in the dark reaction of photosynthesis and is adapted from (Becskei and Serrano, 2000). Every gene has 2 steady state protein concentration predictions with respect to their mRNA expression level in ambient and elevated CO2. The ratio between the steady state protein level at elevated and ambient CO2 was used as one of the inputs for the photosynthesis metabolic model in order to account for changes in gene expression of these metabolic genes in response to CO2 levels. The model assumes constant temperature and constant concentration of RNA polymerase (El Samad et al, 2005). 

## References

Becskei, A., & Serrano, L. (2000). Engineering stability in gene networks by autoregulation. Nature, 405(6786), 590.

El Samad, H., Khammash, M., Petzold, L., & Gillespie, D. (2005). Stochastic modelling of gene regulatory networks. International Journal of Robust and Nonlinear Control, 15(15), 691-711.

## Running the model

```
python src/GrMC.py
```

## Model Inputs/Outputs

### Inputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
Initial protein content | g/L | float | Initial protein concentration before the model is executed.
L | per day | float | Protein synthesis rate.
U | per day | float | Protein degradation rate.
d | NA | float | Effectiveness factor of feedback loop into protein translation. 
r | NA | float | Normalized mRNA expression level from microarray data.

### Outputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
p | g/L | float | Predicted steady state protein level.


