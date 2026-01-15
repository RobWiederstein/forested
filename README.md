# Optimism Bias in Ecological Modeling

## Project Overview

This project explores the hazards of ignoring **spatial autocorrelation** in ecological modeling. Using the `forested` package and forest structure data from Washington State, this Quarto presentation demonstrates how standard random cross-validation yields overly optimistic performance estimates by allowing models to "cheat" via nearby neighbors. The analysis utilizes the `spatialsample` package to visualize and compare three distinct validation strategies—**Random** (the baseline), **Spatial Blocking** (geographic separation), and **Environmental Clustering** (ecological separation)—to establish robust, geographically transferable model performance metrics.