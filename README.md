# Coffee-SDM-Analysis: The R script follows the following steps to perform species distribution modelling for Coffee arabica.
## Load required R packages: 
loading the necessary R libraries that are essential for handling spatial data, running models, and evaluating results.
## Load study area shapefile: 
loading the study area shapefile to later mask all predictor variables to the extent of the study area
## Import species occurrence data: 
Species presence records (Arabica coffee occurrence location GPS records) are read into R from a CSV file in the computer drive. The CSV file contains geographic coordinates (latitude and longitude) where the species has been observed.
## Load and prepare environmental predictor variables:
The 19 Bioclimatic variables (temperature, precipitation-related) and other environmental variables (soil, topographic) are loaded as raster layers and stacked into a single object.
## Data preprocessing:
Spatial filtering to remove duplicate points, checks for multicollinearity among predictor variables using the variance inflation factor (VIF), and selects the predictor variables with no collinearity problem.
## Partition data for model training and testing:
Occurrence data are split into training and testing subsets to enable robust model evaluation.
## Model fitting and prediction: 
Species distribution models (MaxEnt, Random Forest, SVM) are fitted using the training data and environmental predictors, and predictions are made to generate suitability maps for the study area. Visualize model performance using shiny app, visualize plots using ROC-AUC (area under curve plots), and plot variable importance values
## Model evaluation and output: 
Model performance is assessed using metrics such as AUC and TSS. The final habitat suitability maps are saved as raster files for visualization and further climate scenarios 
## Area changes statistics analysis:
Predicted habitat suitability maps are converted into binary maps (suitable/unsuitable) using optimal thresholds. The script then calculates the total suitable area for the current and each future scenario, and determines the niche change (gain, loss, or remaining stable).
## Output and visualization: 
The variable importance plots, model evaluation results, final suitability maps, area change statistics, and niche change plots are saved for reporting and publication under each step.
## EUDR Comparison:
In this step, the projected coffee suitability maps are compared with the European Union Deforestation Regulation (EUDR) compliance requirements. The analysis overlays the projected coffee suitability maps with the global forest cover map used by EUDR to identify areas that are both climatically suitable and meet deforestation-free criteria, to assess how much of the future suitable area aligns with EUDR compliance.
## Remark
## outputs uploaded to the GitHub repository
-Projected suitability maps for five study sites, two scenarios(SSSP245 and SSP585), four time periods(2023s,2050s,2070, and 2090s)
## processed soil data
The final processed soil data can be provided upon request due to the large file size.
### Coffee occurrence locations
Coffee occurrence locations are uploaded here in shapefile format to visualize in GIS software.
