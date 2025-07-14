## Overview
This repository contains an ecological statistical model for determining optimal forest ranger staffing levels in protected areas. The model integrates landscape ecology, threat assessment, and conservation priorities to calculate ranger requirements that balance ecological protection with operational efficiency.

## Key Features
1. **Multi-scale Analysis**: Incorporates both landscape-level and site-specific predictors
2. **Spatial Modeling**: Uses landscape metrics to quantify habitat patterns
3. **Uncertainty Quantification**: Provides confidence intervals for all predictions
4. **Adaptive Framework**: Model parameters can be updated with new field data
5. **Visual Reporting**: Generates comprehensive visualizations for stakeholders

## Methodology
The model uses a hierarchical Bayesian approach with the following structure:

## Requirements
- R (v4.1.0+)
- `raster`, `sf` (spatial analysis)
- `glmmTMB`, `MuMIn` (statistical modeling)
- `ggplot2`, `patchwork` (visualization)
- `DHARMa` (model diagnostics)
