# Protein translation model (PTM)

## Model Description

This model essentially predicts the steady state protein concentration of selected enzymes in the dark reaction of photosynthesis and is adapted from (Becskei and Serrano, 2000). Every gene has 2 steady state protein concentration predictions with respect to their mRNA expression level in ambient and elevated CO2. The ratio between the steady state protein level at elevated and ambient CO2 was used as one of the inputs for the photosynthesis metabolic model in order to account for changes in gene expression of these metabolic genes in response to CO2 levels. The model assumes constant temperature and constant concentration of RNA polymerase (El Samad et al, 2005). 

## References

Becskei, A., & Serrano, L. (2000). Engineering stability in gene networks by autoregulation. Nature, 405(6786), 590.

El Samad, H., Khammash, M., Petzold, L., & Gillespie, D. (2005). Stochastic modelling of gene regulatory networks. International Journal of Robust and Nonlinear Control, 15(15), 691-711.

## Running the model

### Without yggdrasil
```
python src/GrCM_plain.py
```

### With yggdrasil
```
yggrun GrCM.yml GrCM_tofile.yml
```

## Model Inputs/Outputs

### Inputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
Initial protein content | g/L | float | Initial protein concentration before the model is executed.
L | per day | float | Protein synthesis rate.
U | per day | float | Protein degradation rate.
d | NA | float | Effectiveness factor of feedback loop into protein translation. 
r | NA | float | Condition specific transcript levels of the gene.

### Outputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
p | g/L | float | Predicted steady state protein level.

## Additional notes regarding reproducibility of some figures and tables pertaining to PTM model in Kannan et al., 2019 (doi:'https://doi.org/10.1093/insilicoplants/diz008')

Supplemental figure S2 provides an example simulation of selection of 'D' factor for two Glyma IDs that is used as one of the input parameters in 'Input/GrCM_input.txt' file. Data pertaining to this figure is provided as an OriginPro (https://www.originlab.com/) data analysis and graphing software file 'SuppFigureS2.opju'

Supplemental table S2 provides some information contained in 'Input/GrCM_input.txt' file

Supplemental table S5 provides information contained in the 'Output/GrCM_output.txt' file
