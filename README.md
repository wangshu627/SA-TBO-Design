# SA-TBO: Space-Filling Design Expansion via Block Design

This repository contains the R implementation of the **SA-TBO (Simulated Annealing with Targeted Bottleneck Optimization)** framework. This methodology is designed for generating and expanding high-dimensional space-filling designs efficiently via block design techniques.

## Overview

Generating large-scale space-filling designs (like Latin Hypercube Designs) can be computationally expensive. This repository provides a two-stage construction framework that:
1. Implements a fast, block-based expansion method (Algorithm 1).
2. Optimizes the block vectors and refines the design using the SA-TBO algorithm to maximize spatial uniformity (measured by $L_1$ and $L_2$ distances).

## Repository Structure

The code is organized into three main analytical scripts:

* `01_core_algorithm_and_scalability.R`
    * **Purpose:** Contains the core Algorithm 1, distance metrics ($L_1$/$L_2$ upper bounds, column correlation), and the implementation of the SA-TBO block vector optimization. Includes batch execution for large-scale scalability testing.
* `02_baseline_comparisons.R`
    * **Purpose:** Compares the proposed SA-TBO framework against established external baselines, including Morris & Mitchell (maximin LHS) and Maximin Sliced LHD (SLHD).
* `03_case_study_borehole.R`
    * **Purpose:** A practical engineering case study applying the expanded designs to the surrogate modeling of the Borehole function. It includes 30 independent trials to ensure statistical significance, calculates RMSE/MAE using Kriging models, and generates high-quality PDF plots for visualization.

## Dependencies

To run the scripts, you will need **R** installed on your system along with the following packages:

```R
install.packages(c("lhs", "SLHD", "DiceKriging"))
```

* `lhs`: Used for generating Morris & Mitchell baseline designs.
* `SLHD`: Used for generating Maximin Sliced Latin Hypercube baseline designs.
* `DiceKriging`: Used for building the surrogate/Kriging models in the Borehole case study.

## Usage

You can run each script independently in your preferred R environment (e.g., RStudio). 

1. **For basic methodology and scalability:** Run `01_core_algorithm_and_scalability.R`. Results will be printed to the console and saved to `scalability_results.csv`.
2. **For comparative analysis:** Run `02_baseline_comparisons.R`. The benchmark results will be exported to `Expanded_Design_Comparisons.csv`.
3. **For the surrogate modeling application:** Run `03_case_study_borehole.R`. This will generate several PDF figures (e.g., `fig_boxplot_rmse5.pdf`, `fig_2d_projection5.pdf`) and a summary CSV.

## Citation

If you find this code or the SA-TBO framework useful in your research, please consider citing our paper (currently under review for *IISE Transactions*):

```bibtex
@article{Wang2026SATBO,
  title={Analysis and Improvement of Space-Filling Design Expansion via Block Design},
  author={Wang, Shu and [Co-authors]},
  journal={IISE Transactions},
  year={2026},
  note={Under review}
}
```
*(Note: Please update the citation details once the paper is published.)*

## License

This project is licensed under the MIT License - see the LICENSE file for details.
