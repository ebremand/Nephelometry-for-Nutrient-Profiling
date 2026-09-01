# Nephelometry for Nutrient Profiling

This repository contains the R scripts used for the data analysis and figure generation associated with the manuscript:
**"Comparative Evaluation of Growth Assessment Methods in Filamentous Fungi Reveals the Potential of Nephelometry for Nutrient Profiling"**

The study evaluates nephelometry as a high-throughput method for monitoring fungal growth and characterizing nutrient utilization in filamentous fungi. Nephelometric measurements are compared with spectrophotometry, colony diameter, and fungal biomass measurements obtained from Petri dishes.

Corresponding author : Justine Colou, justine.colou@univ-angers.fr

## Repository contents

### `Growth_Curve_Analysis.R`

This script processes nephelometric (NTU) and spectrophotometric (OD) growth kinetics and calculates the growth parameters used in the study.

Raw growth curves were subjected to three preprocessing steps prior to parameter extraction. First, values recorded before the minimum signal detected for each growth curve were replaced by this minimum value to correct potential condensation-related artifacts occurring at the beginning of measurements. Second, baseline correction was performed by subtracting the minimum value from each curve, setting the initial signal to zero. Third, corrected curves were smoothed using a LOESS (LOcally Estimated Scatterplot Smoothing) regression to reduce experimental noise while preserving the overall growth dynamics.

The growth rate was estimated from the highest slope observed along each growth curve. A sliding window approach was applied using a 10 h window over the growth phase preceding the maximum signal intensity. For each window, a linear regression based on the least-squares method was performed to estimate the slope of signal increase over time. The highest slope among all tested windows was considered the maximum growth rate and was expressed as NTU·h⁻¹ or OD·h⁻¹.

The regression line corresponding to the growth rate was then used to estimate additional kinetic parameters (Figure 1). The lag time was calculated as the time point at which the tangent line reached the baseline, corresponding to the moment when the extrapolated growth signal reached zero, and was expressed in hours. The maximum carrying capacity (MCC) was defined as the maximum signal value experimentally reached by each growth curve and was expressed as NTU or OD. The MCC time was estimated as the intersection between the tangent line corresponding to the growth rate and the maximum experimentally observed signal, and was expressed in hours. The exponential duration was calculated as the time interval between the lag time and the MCC time, and was expressed in hours. Finally, the area under the growth curve (AUC) was calculated by numerical integration using the trapezoidal method implemented in the trapz function from the pracma R package v2.4.6 and was expressed as NTU·h or OD·h.


### `Figure_Publication.R`

This script generates the figures presented in the manuscript and supplementary material, including:

* Spearman correlation matrices between growth parameters and dry weight
* Comparisons between nephelometry, spectrophotometry, colony diameter, and dry biomass
* Standardized boxplots of growth measurements
* Coefficient of variation analyses
* Growth curves obtained by nephelometry and spectrophotometry
* Correlation analyses after extended growth periods
* Evaluation of nutrient concentration conditions
* Correlation analyses following nutrient optimization

## Data

The datasets required to run the analyses are provided as Supplementary Material of the associated publication.
- Supplementary Dataset 1: Growth measurements of A. brassicicola on Petri dishes on 12 nutrients. Colony diameter is expressed in mm, and fresh and dry weight are expressed in mg. Data used to generate Figures 2, 3, 4, 5, 6, and 8.
- Supplementary Dataset 2: Growth kinetic values of A. brassicicola on 12 nutrients measured by spectrophotometry (OD) and nephelometry (NTU). Data used to generate Figures 2, 3, 4, 5, and 6.
- Supplementary Dataset 3: Growth kinetic values of A. brassicicola measured by nephelometry (NTU) at different glucose concentrations (in g·L-1) and asparagine concentrations (in mM). Data used to generate Figure 7.
- Supplementary Dataset 4: Growth kinetic values of A. brassicicola measured by nephelometry (NTU) on 12 nutrients using optimized nitrogen and glucose concentrations (Figure 7). Data used to generate Figure 8.
