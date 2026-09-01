# Nephelometry for Nutrient Profiling

This repository contains the R scripts used for the data analysis and figure generation associated with the manuscript:
**"Comparative Evaluation of Growth Assessment Methods in Filamentous Fungi Reveals the Potential of Nephelometry for Nutrient Profiling"**

The study evaluates nephelometry as a high-throughput method for monitoring fungal growth and characterizing nutrient utilization in filamentous fungi. Nephelometric measurements are compared with spectrophotometry, colony diameter, and fungal biomass measurements obtained from Petri dishes.

## Repository contents

### `Nephelometry_Analysis.R`

This script contains the analyses performed on the nephelometric data, including data processing and growth parameter calculations.

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
