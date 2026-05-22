# Tombe & Zhu (2019) Replication

## Authors

Zhuangfu Zhao & Tianyi Liang

## Original Paper

Tombe, Trevor, and Xiaodong Zhu. 2019. "Trade, Migration, and Productivity: A Quantitative Analysis of China." *American Economic Review* 109 (5): 1843–72.

## About This Replication

The original replication package was written in MATLAB. In this repository, we:

1. **Trace the theoretical foundations** of Tombe & Zhu (2019): The model builds on the quantitative trade framework of Eaton & Kortum (2002), extended by Caliendo & Parro (2015) to incorporate input-output linkages. Tombe & Zhu (2019) further introduce internal migration and land usage into this framework.
2. **Replicate the model equilibrium computation in Python**.
3. **Provide detailed annotations** that explain the code logic and its correspondence to the paper.

## Repository Structure

- `0_Fréchet/` — Basic properties of the Fréchet and Gumbel distributions.
- `1_Eaton-Kortum/` — Model derivation of Eaton & Kortum (2002) and code examples.
- `2_Caliendo-Parro/` — Model derivation of Caliendo & Parro (2015) and code examples.
- `3_Tombe-Zhu/` — Model derivation of Tombe & Zhu (2019) and Python replication of the original MATLAB code.
- `origin_matlab/` — Original replication code from Tombe & Zhu (MATLAB/Stata).
- `data/` — Data files used in the replication.

Within each folder, files named `main_*` are the primary notebooks; `module_*` are the modules called by the `main_*` notebooks; `module_XX_interpretation` files are explanatory walkthroughs for the corresponding modules, provided as notebooks for easy cross-referencing between code and model derivation, and should not be executed.

## References

&nbsp;&nbsp;&nbsp;&nbsp;Caliendo, Lorenzo, and Fernando Parro. 2015. "Estimates of the Trade and Welfare Effects of NAFTA." *Review of Economic Studies* 82 (1): 1–44.

&nbsp;&nbsp;&nbsp;&nbsp;Dekle, Robert, Jonathan Eaton, and Samuel Kortum. 2008. "Global Rebalancing with Gravity: Measuring the Burden of Adjustment." *IMF Staff Papers* 55 (3): 511–40.

&nbsp;&nbsp;&nbsp;&nbsp;Eaton, Jonathan, and Samuel Kortum. 2002. "Technology, Geography, and Trade." *Econometrica* 70 (5): 1741–79.

&nbsp;&nbsp;&nbsp;&nbsp;Kortum, Samuel S. 1997. "Research, Patenting, and Technological Change." *Econometrica* 65 (6): 1389.
