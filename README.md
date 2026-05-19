# Tombe & Zhu (2019) Replication

## Authors

Zhuangfu Zhao & Tianyi Liang

## Original Paper

Tombe, Trevor, and Xiaodong Zhu. 2019. "Trade, Migration, and Productivity: A Quantitative Analysis of China." *American Economic Review* 109 (5): 1843–72.

## About This Replication

The original replication package was written in MATLAB. In this repository, we:

1. **Trace the theoretical foundations** of Tombe & Zhu (2019): The model builds on the quantitative trade framework of Eaton & Kortum (2002), extended by Caliendo & Parro (2015) to incorporate input-output linkages. Tombe & Zhu (2019) further introduce internal migration and land use into this framework.
2. **Replicate the model equilibrium computation in Python**.
3. **Provide detailed annotations** that explain the code logic and its correspondence to the equations in the paper.

## Repository Structure

- `origin_matlab/` — Original replication code from Tombe & Zhu (MATLAB/Stata).
- `data/` — Data files used in the replication.
- `TombeZhu2019.ipynb` — Main notebook covering: the theoretical foundations (EK and CP models), the key content of Tombe & Zhu (2019), and the Python replication with detailed annotations.
- `main_simulate.py` — Contains the `eqm()` function called by `TombeZhu2019.ipynb`, encoding the model's equilibrium conditions (the system of equations to be solved).
- `main_simulate_interpretation.ipynb` — Line-by-line interpretation of `main_simulate.py`; provided as a notebook for readability but should not be executed.
- `Frechet.py` — Functions called by `TombeZhu2019.ipynb` to illustrate basic properties of the Fréchet and Gumbel distributions.

## References

&nbsp;&nbsp;&nbsp;&nbsp;Eaton, Jonathan, and Samuel Kortum. 2002. "Technology, Geography, and Trade." *Econometrica* 70 (5): 1741–79.

&nbsp;&nbsp;&nbsp;&nbsp;Dekle, Robert, Jonathan Eaton, and Samuel Kortum. 2008. "Global Rebalancing with Gravity: Measuring the Burden of Adjustment." *IMF Staff Papers* 55 (3): 511–40.

&nbsp;&nbsp;&nbsp;&nbsp;Caliendo, Lorenzo, and Fernando Parro. 2015. "Estimates of the Trade and Welfare Effects of NAFTA." *Review of Economic Studies* 82 (1): 1–44.
