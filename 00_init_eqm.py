"""
00_init_eqm.py
==============
This script sets up the initial equilibrium of the model from Tombe & Zhu (2019)
"Trade, Migration and Productivity: A Quantitative Analysis of China". 

Data is loaded from data/raw/.
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import jax
import jax.numpy as jnp

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(SCRIPT_DIR, "data", "raw")


# ===========================================================================
# 1. Merge in the mu and price data
# ===========================================================================
data_emp = pd.read_csv(os.path.join(RAW_DIR, "employment_realGDP_data.csv"))
region = data_emp["province"].tolist()
N = len(region) + 1  # 30 provinces + 1 international = 31

rYLa_data = data_emp["rYLa2000"].to_numpy().reshape(-1, 1)   # data.data(:,1)
rYLn_data = data_emp["rYLn2000"].to_numpy().reshape(-1, 1)   # data.data(:,2)
La_data   = data_emp["La2000"].to_numpy().reshape(-1, 1)     # data.data(:,3)
Ln_data   = data_emp["Ln2000"].to_numpy().reshape(-1, 1)     # data.data(:,4)
l = (La_data + Ln_data)

# L(1:N,1) = sum(La_data + Ln_data)  -> all entries the same total
L = np.full((N, 1), float((La_data + Ln_data).sum()))

# Real income per worker (province-level mean)
rYL_data = (rYLa_data * La_data + rYLn_data * Ln_data) / (La_data + Ln_data)

# mu = [rYLa_data rYLn_data]' ; mu = mu(:)
# stacks rows interleaved: prov1_ag, prov1_na, prov2_ag, prov2_na, ...
mu = np.column_stack([rYLa_data, rYLn_data]).reshape(-1, 1)

rYLa_data2005 = data_emp["rYLa2005"].to_numpy().reshape(-1, 1)  # data.data(:,5)
rYLn_data2005 = data_emp["rYLn2005"].to_numpy().reshape(-1, 1)  # data.data(:,6)
drYLa_data    = data_emp["drYLa_data"].to_numpy().reshape(-1, 1)  # data.data(:,7)
drYLn_data    = data_emp["drYLn_data"].to_numpy().reshape(-1, 1)  # data.data(:,8)
drYL_data = np.column_stack([drYLa_data, drYLn_data]).reshape(-1, 1)

# International values:  L(N,1) = 2102.978587
L[N - 1, 0] = 2102.978587


# ===========================================================================
# 2. Load Initial Trade Shares
#    (lines 30-33 of initial_eqm.m)
# ===========================================================================
pi_ag = pd.read_csv(os.path.join(RAW_DIR, "trade_ag.csv"), header=None).to_numpy()
pi_ag = pi_ag / pi_ag.sum(axis=1, keepdims=True)
pi_na = pd.read_csv(os.path.join(RAW_DIR, "trade_na.csv"), header=None).to_numpy()
pi_na = pi_na / pi_na.sum(axis=1, keepdims=True)


# ===========================================================================
# 3. Load Changes in Trade Costs
#    (lines 36-40 of initial_eqm.m)
# ===========================================================================
data_tau = pd.read_csv(os.path.join(RAW_DIR, "tauhat.csv"))
# In MATLAB, reshape(data.data(:,1), N, N)' produces an N x N matrix.
# The MATLAB column vector is column-major so reshape(N,N)' converts to
# the equivalent (i,j) interpretation. In NumPy this is reshape((N,N)).
dni_ag_measured = data_tau["dni_ag"].to_numpy().reshape((N, N))
dni_na_measured = data_tau["dni_na"].to_numpy().reshape((N, N))
dni_asym_ag     = data_tau["dni_asym_ag"].to_numpy().reshape((N, N))
dni_asym_na     = data_tau["dni_asym_na"].to_numpy().reshape((N, N))


# ===========================================================================
# 4. Load Migration Shares
#    (lines 44-47 of initial_eqm.m)
# ===========================================================================
data_mij2000 = pd.read_csv(os.path.join(RAW_DIR, "mij2000.csv"))
# reshape(data.data(:,2), 60, 60)'  -> 60x60 matrix, each row = origin
mij2000 = data_mij2000["mij"].to_numpy().reshape((60, 60))

data_mij2005 = pd.read_csv(os.path.join(RAW_DIR, "mij2005.csv"))
mij2005 = data_mij2005["mij"].to_numpy().reshape((60, 60))


# ===========================================================================
# 5. Initial Employment Vector & Hukou Registrations
#    (lines 50-58 of initial_eqm.m)
# ===========================================================================
# L1=[La_data Ln_data]'; L1=L1(:)  -> stack interleaved (only the 30 provinces)
L1 = np.column_stack([La_data, Ln_data]).reshape(-1, 1)  # 60x1

# L0 = (L1' * inv(mij2000))'   - hukou registrations
L0 = np.linalg.solve(mij2000.T, L1).reshape(-1, 1)

# reshape(L1', 2, N-1)'   -> (N-1) x 2
temp = L1.reshape(-1, 2)
La = temp[:, 0:1]
Ln = temp[:, 1:2]


# ===========================================================================
# 6. Production Function Input Shares
#    (lines 61-77 of initial_eqm.m)
# ===========================================================================
beta_ag = 0.287
beta_na = 0.219    # Brandt, Restuccia, Adamopolous + IO data
eta_ag  = 0.277
eta_na  = 0.025

# Intermediate Input Shares (gamma matrix)
gamma = np.array([[0.397, 0.603],
                  [0.063, 0.937]])

# Common Parameters
alpha   = 1 - 0.13     # housing share => 0.87
kappa   = 1.5
theta   = 4.0
epsilon = 0.09545      # agriculture's share of total end-use


# ===========================================================================
# 7. Initial total nominal expenditures and income (fixed-point loop)
#    (lines 87-115 of initial_eqm.m)
# ===========================================================================
I0    = np.zeros((N, 1))
I     = np.ones((N, 1))
X0_ag = np.zeros((N, 1))
X_ag  = np.ones((N, 1))
X0_na = np.zeros((N, 1))
X_na  = np.ones((N, 1))

# Pre-declare variables from inside the loop so they exist after it
R_ag = R_na = Ir = Iu = None

while ((X0_ag - X_ag) ** 2).sum() + ((X0_na - X_na) ** 2).sum() > 1e-20:
    X0_ag = X_ag.copy()
    X0_na = X_na.copy()

    R_ag = pi_ag.T @ X0_ag
    R_na = pi_na.T @ X0_na

    Ir = (beta_ag + eta_ag) * R_ag / alpha   # Rural income
    Iu = (beta_na + eta_na) * R_na / alpha   # Urban income
    I = Ir + Iu                              # Total provincial income
    Ir = Ir / I.sum()
    Iu = Iu / I.sum()

    Da_r = alpha * epsilon * Ir
    Da_u = alpha * epsilon * Iu
    Dn_r = alpha * (1 - epsilon) * Ir
    Dn_u = alpha * (1 - epsilon) * Iu
    Ds_r = (1 - alpha) * Ir
    Ds_u = (1 - alpha) * Iu

    X_ag = (Da_r + Da_u
            + (1 - beta_ag - eta_ag) * gamma[0, 0] * R_ag
            + (1 - beta_na - eta_na) * gamma[1, 0] * R_na)
    X_na = (Dn_r + Dn_u
            + (1 - beta_na - eta_na) * gamma[1, 1] * R_na
            + (1 - beta_ag - eta_ag) * gamma[0, 1] * R_ag)


# ===========================================================================
# 8. Initial and 2005 real income per worker
#    (lines 120-131 of initial_eqm.m)
# ===========================================================================
# Vi has 2*(N-1) rows  (only the 30 provinces, two sectors)
Via = rYLa_data.copy()
Vin = rYLn_data.copy()
Via2005 = rYLa_data * drYLa_data
Vin2005 = rYLn_data * drYLn_data

Vi      = np.column_stack([Via,       Vin]).reshape(-1, 1)
Vi2005  = np.column_stack([Via2005, Vin2005]).reshape(-1, 1)

# dVi_data = [Vi2005./Vi; 1.055; 1.055]
dVi_data = np.vstack([Vi2005 / Vi, np.array([[1.055], [1.055]])])
temp = dVi_data.reshape(-1, 2)
dVa_data = temp[:, 0:1]
dVn_data = temp[:, 1:2]


# ===========================================================================
# 9. Initial migration costs (per effective labour)
#    (lines 134-145 of initial_eqm.m)
# ===========================================================================
TWO_NM1 = 2 * (N - 1)  # 60

Vimat = np.tile(Vi, (1, TWO_NM1))
Vjmat = np.tile(Vi.T, (TWO_NM1, 1))

cij = (mij2000 / np.tile(np.diag(mij2000).reshape(-1, 1), (1, TWO_NM1))) ** (1 / kappa) \
      / (Vjmat / Vimat)
cij2000 = cij.copy()

Vimat2005 = np.tile(Vi2005, (1, TWO_NM1))
Vjmat2005 = np.tile(Vi2005.T, (TWO_NM1, 1))
cij2005 = (mij2005 / np.tile(np.diag(mij2005).reshape(-1, 1), (1, TWO_NM1))) ** (1 / kappa) \
          / (Vjmat2005 / Vimat2005)

# Recompute mij from cij2000 and Vjmat
mij_local = (cij2000 * Vjmat) ** kappa
mij_local = mij_local / mij_local.sum(axis=1, keepdims=True)
cijhat = np.ones((TWO_NM1, TWO_NM1))
L1base = (L0.T @ mij_local).T            # 60x1
L2005  = (L0.T @ mij2005).T              # 60x1


# ===========================================================================
# 10. Land Rebate Adjustments to Migration Costs 1111
#     (lines 149-169 of initial_eqm.m)
# ===========================================================================
M_mL_term = L1base / (np.diag(mij_local).reshape(-1, 1) * L0)
temp = M_mL_term.reshape(-1, 2)
M_mL_term_ag = temp[:, 0:1]
M_mL_term_na = temp[:, 1:2]

Cnnjj_ag = 1 + (eta_ag + (1 - alpha) * beta_ag) * M_mL_term_ag / (alpha * beta_ag)
Cnnjj_na = 1 + (eta_na + (1 - alpha) * beta_na) * M_mL_term_na / (alpha * beta_na)
Cnnjj = np.column_stack([Cnnjj_ag, Cnnjj_na]).reshape(-1, 1)

Cnijk = (1 - np.eye(TWO_NM1)) * np.tile(Cnnjj, (1, TWO_NM1)) * cij
Cnijk2000 = Cnijk.copy()
Cnnjj2000 = Cnnjj.copy()

# 2005
M_mL_term2005 = L2005 / (np.diag(mij2005).reshape(-1, 1) * L0)
temp = M_mL_term2005.reshape(-1, 2)
M_mL_term_ag2005 = temp[:, 0:1]
M_mL_term_na2005 = temp[:, 1:2]
Cnnjj_ag2005 = 1 + (eta_ag + (1 - alpha) * beta_ag) * M_mL_term_ag2005 / (alpha * beta_ag)
Cnnjj_na2005 = 1 + (eta_na + (1 - alpha) * beta_na) * M_mL_term_na2005 / (alpha * beta_na)
Cnnjj2005 = np.column_stack([Cnnjj_ag2005, Cnnjj_na2005]).reshape(-1, 1)
Cnijk2005 = (1 - np.eye(TWO_NM1)) * np.tile(Cnnjj2005, (1, TWO_NM1)) * cij2005


# ===========================================================================
# 11. International initial labour distribution
#     (lines 172-173 of initial_eqm.m)
# ===========================================================================
# Ln_data and La_data here only have N-1 rows. Append the international row.
Ln_data_full = np.vstack([Ln_data, np.zeros((1, 1))])
La_data_full = np.vstack([La_data, np.zeros((1, 1))])

Ln_data_full[N - 1, 0] = L[N - 1, 0] / (1 + (beta_ag / beta_na) * (R_ag[N - 1, 0] / R_na[N - 1, 0]))
La_data_full[N - 1, 0] = L[N - 1, 0] - Ln_data_full[N - 1, 0]


# ===========================================================================
# 12. Average symmetric trade costs
#     (lines 176-177 of initial_eqm.m)
# ===========================================================================
diag_pi_ag = np.diag(pi_ag).reshape(-1, 1)
diag_pi_na = np.diag(pi_na).reshape(-1, 1)
tau_ag = ((diag_pi_ag @ diag_pi_ag.T) / (pi_ag * pi_ag.T)) ** (1.0 / (2.0 * theta))
tau_na = ((diag_pi_na @ diag_pi_na.T) / (pi_na * pi_na.T)) ** (1.0 / (2.0 * theta))


# ===========================================================================
# 13. Get the Initial Equilibrium of the Model (solve the system)
#     (lines 181-191 of initial_eqm.m)
# ===========================================================================
dTa = np.ones((N, 1))
dTn = np.ones((N, 1))
dni_ag = np.ones((N, N))
dni_na = np.ones((N, N))


def main_simulate(X, params):
    """
    Python port of main_simulate.m
    Returns the residual vector F that fsolve drives to zero.
    """
    Np         = params["N"]
    theta_p    = params["theta"]
    beta_ag_p  = params["beta_ag"]
    beta_na_p  = params["beta_na"]
    eta_ag_p   = params["eta_ag"]
    eta_na_p   = params["eta_na"]
    dni_ag_p   = params["dni_ag"]
    dni_na_p   = params["dni_na"]
    pi_ag_p    = params["pi_ag"]
    pi_na_p    = params["pi_na"]
    L0_p       = params["L0"]
    kappa_p    = params["kappa"]
    alpha_p    = params["alpha"]
    dTa_p      = params["dTa"]
    dTn_p      = params["dTn"]
    Vi_p       = params["Vi"]
    mij2000_p  = params["mij2000"]
    R_ag_p     = params["R_ag"]
    R_na_p     = params["R_na"]
    Ir_p       = params["Ir"]
    Iu_p       = params["Iu"]
    epsilon_p  = params["epsilon"]
    La_data_p  = params["La_data"]
    Ln_data_p  = params["Ln_data"]
    L1base_p   = params["L1base"]
    gamma_p    = params["gamma"]
    Cnijk_p    = params["Cnijk"]

    X = X.reshape(-1, 1)
    dwa = X[0:Np]
    dwn = X[Np:2 * Np]
    dPa = X[2 * Np:3 * Np]
    dPn = X[3 * Np:4 * Np]
    dLa = X[4 * Np:5 * Np]
    dLn = X[5 * Np:6 * Np]
    dL = np.column_stack([dLa, dLn]).reshape(-1, 1)

    # New sectoral revenues
    R_ag_new = dwa * dLa * R_ag_p
    R_na_new = dwn * dLn * R_na_p

    # New household final demand
    Ir_new = (beta_ag_p + eta_ag_p) * R_ag_new / alpha_p
    dIr = Ir_new / Ir_p
    Iu_new = (beta_na_p + eta_na_p) * R_na_new / alpha_p
    dIu = Iu_new / Iu_p
    Da_new_r = alpha_p * epsilon_p * Ir_new
    Da_new_u = alpha_p * epsilon_p * Iu_new
    Dn_new_r = (alpha_p * (1 - epsilon_p)) * Ir_new
    Dn_new_u = (alpha_p * (1 - epsilon_p)) * Iu_new

    # Change in land prices
    dra = dwa * dLa
    drn = dwn * dLn

    # New trade shares
    inside_ag = (dwa.T ** beta_ag_p
                 * (dPa.T ** gamma_p[0, 0] * dPn.T ** gamma_p[0, 1]) ** (1 - beta_ag_p - eta_ag_p)
                 * dra.T ** eta_ag_p
                 / (dTa_p.T ** (1.0 / theta_p)))
    kron_inside_ag = np.kron(inside_ag, np.ones((Np, 1)))
    denominator_ag = np.sum(pi_ag_p * (dni_ag_p * kron_inside_ag) ** (-theta_p), axis=1, keepdims=True)
    pi_ag_new = (pi_ag_p * (dni_ag_p * kron_inside_ag) ** (-theta_p)) / np.kron(denominator_ag, np.ones((1, Np)))

    inside_na = (dwn.T ** beta_na_p
                 * (dPa.T ** gamma_p[1, 0] * dPn.T ** gamma_p[1, 1]) ** (1 - beta_na_p - eta_na_p)
                 * drn.T ** eta_na_p
                 / (dTn_p.T ** (1.0 / theta_p)))
    kron_inside_na = np.kron(inside_na, np.ones((Np, 1)))
    denominator_na = np.sum(pi_na_p * (dni_na_p * kron_inside_na) ** (-theta_p), axis=1, keepdims=True)
    pi_na_new = (pi_na_p * (dni_na_p * kron_inside_na) ** (-theta_p)) / np.kron(denominator_na, np.ones((1, Np)))

    # Sector expenditures
    X_ag_new = np.linalg.solve(pi_ag_new.T, R_ag_new)
    X_na_new = np.linalg.solve(pi_na_new.T, R_na_new)

    # Real income per effective labour
    dP = dPa ** epsilon_p * dPn ** (1 - epsilon_p)
    dVa = dwa / (dP ** alpha_p * dra ** (1 - alpha_p))
    dVn = dwn / (dP ** alpha_p * drn ** (1 - alpha_p))
    dV = np.column_stack([dVa, dVn]).reshape(-1, 1)

    # Migration shares (inner fixed-point loop)
    two_nm1 = 2 * (Np - 1)
    mnn1 = np.diag(mij2000_p).reshape(-1, 1)
    mnn0 = np.ones((two_nm1, 1))
    mij_loc = None
    while ((mnn0 - mnn1) ** 2).sum() > 1e-20:
        mnn0 = mnn1.copy()
        M_mL_term = (L1base_p * dL[0:two_nm1]) / (mnn0 * L0_p)
        temp_local = M_mL_term.reshape(-1, 2)
        M_mL_term_ag_l = temp_local[:, 0:1]
        M_mL_term_na_l = temp_local[:, 1:2]
        Cnnjj_ag_l = 1 + (eta_ag_p + (1 - alpha_p) * beta_ag_p) * M_mL_term_ag_l / (alpha_p * beta_ag_p)
        Cnnjj_na_l = 1 + (eta_na_p + (1 - alpha_p) * beta_na_p) * M_mL_term_na_l / (alpha_p * beta_na_p)
        Cnnjj_l = np.column_stack([Cnnjj_ag_l, Cnnjj_na_l]).reshape(-1, 1)
        cij_l = Cnijk_p + np.eye(two_nm1) * np.tile(Cnnjj_l, (1, two_nm1))
        Vjmat_l = np.tile((dV[0:two_nm1] * Vi_p).T, (two_nm1, 1))
        mij_loc = (cij_l * Vjmat_l) ** kappa_p
        mij_loc = mij_loc / mij_loc.sum(axis=1, keepdims=True)
        mnn1 = np.mean(np.column_stack([np.diag(mij_loc).reshape(-1, 1), mnn0]), axis=1, keepdims=True)

    # New labour allocation
    L1_loc = (L0_p.T @ mij_loc).T

    # System of equations
    F = np.vstack([
        dPa - denominator_ag ** (-1.0 / theta_p),
        dPn - denominator_na ** (-1.0 / theta_p),
        L1_loc[0:two_nm1] / L1base_p[0:two_nm1] - dL[0:two_nm1],
        np.array([[(dLa[Np - 1, 0] * La_data_p[Np - 1, 0] + dLn[Np - 1, 0] * Ln_data_p[Np - 1, 0]
                    - (La_data_p[Np - 1, 0] + Ln_data_p[Np - 1, 0]))]]),
        np.array([[dwa[Np - 1, 0] - dwn[Np - 1, 0]]]),
        X_ag_new - (Da_new_r + Da_new_u
                    + (1 - beta_ag_p - eta_ag_p) * gamma_p[0, 0] * R_ag_new
                    + (1 - beta_na_p - eta_na_p) * gamma_p[1, 0] * R_na_new),
        X_na_new - (Dn_new_r + Dn_new_u
                    + (1 - beta_na_p - eta_na_p) * gamma_p[1, 1] * R_na_new
                    + (1 - beta_ag_p - eta_ag_p) * gamma_p[0, 1] * R_ag_new),
    ])
    return F.flatten()


params = {
    "N": N, "theta": theta,
    "beta_ag": beta_ag, "beta_na": beta_na,
    "eta_ag": eta_ag, "eta_na": eta_na,
    "dni_ag": dni_ag, "dni_na": dni_na,
    "pi_ag": pi_ag, "pi_na": pi_na,
    "L0": L0, "kappa": kappa, "alpha": alpha,
    "dTa": dTa, "dTn": dTn,
    "Vi": Vi, "mij2000": mij2000,
    "R_ag": R_ag, "R_na": R_na, "Ir": Ir, "Iu": Iu,
    "epsilon": epsilon,
    "La_data": La_data_full, "Ln_data": Ln_data_full,
    "L1base": L1base, "gamma": gamma, "Cnijk": Cnijk,
}

x0 = np.ones(6 * N)
new = fsolve(lambda x: main_simulate(x, params), x0, full_output=False, xtol=1e-10)

dwa = new[0:N].reshape(-1, 1)
dwn = new[N:2 * N].reshape(-1, 1)
dPa = new[2 * N:3 * N].reshape(-1, 1)
dPn = new[3 * N:4 * N].reshape(-1, 1)
dLa = new[4 * N:5 * N].reshape(-1, 1)
dLn = new[5 * N:6 * N].reshape(-1, 1)
dL = np.column_stack([dLa, dLn]).reshape(-1, 1)
dw = np.column_stack([dwa, dwn]).reshape(-1, 1)


# ===========================================================================
# 14. Post-solve quantities (lines 196-233 of initial_eqm.m)
# ===========================================================================
R_ag_new = dwa * dLa * R_ag
R_na_new = dwn * dLn * R_na

Ir_new = (beta_ag + eta_ag) * R_ag_new / alpha
dIr = Ir_new / Ir
Iu_new = (beta_na + eta_na) * R_na_new / alpha
dIu = Iu_new / Iu
Da_new_r = alpha * epsilon * Ir_new
Da_new_u = alpha * epsilon * Iu_new
Dn_new_r = (alpha * (1 - epsilon)) * Ir_new
Dn_new_u = (alpha * (1 - epsilon)) * Iu_new
Ds_new_r = (1 - alpha) * Ir_new
Ds_new_u = (1 - alpha) * Iu_new

dra = dwa * dLa
drn = dwn * dLn

inside_ag = (dwa.T ** beta_ag
             * (dPa.T ** gamma[0, 0] * dPn.T ** gamma[0, 1]) ** (1 - beta_ag - eta_ag)
             * dra.T ** eta_ag
             / (dTa.T ** (1.0 / theta)))
kron_inside_ag = np.kron(inside_ag, np.ones((N, 1)))
denominator_ag = np.sum(pi_ag * (dni_ag * kron_inside_ag) ** (-theta), axis=1, keepdims=True)
pi_ag_new = (pi_ag * (dni_ag * kron_inside_ag) ** (-theta)) / np.kron(denominator_ag, np.ones((1, N)))

inside_na = (dwn.T ** beta_na
             * (dPa.T ** gamma[1, 0] * dPn.T ** gamma[1, 1]) ** (1 - beta_na - eta_na)
             * drn.T ** eta_na
             / (dTn.T ** (1.0 / theta)))
kron_inside_na = np.kron(inside_na, np.ones((N, 1)))
denominator_na = np.sum(pi_na * (dni_na * kron_inside_na) ** (-theta), axis=1, keepdims=True)
pi_na_new = (pi_na * (dni_na * kron_inside_na) ** (-theta)) / np.kron(denominator_na, np.ones((1, N)))

X_ag_new = np.linalg.solve(pi_ag_new.T, R_ag_new)
X_na_new = np.linalg.solve(pi_na_new.T, R_na_new)

dP = dPa ** epsilon * dPn ** (1 - epsilon)
dVa = dwa / (dP ** alpha * dra ** (1 - alpha))
dVn = dwn / (dP ** alpha * drn ** (1 - alpha))
dV = np.column_stack([dVa, dVn]).reshape(-1, 1)
dPa_vector = np.column_stack([dPa, dPa]).reshape(-1, 1)
dv_vector = np.column_stack([dIr / dLa, dIu / dLn]).reshape(-1, 1)


# ===========================================================================
# 15. Migration Shares (post-equilibrium fixed-point) lines 236-251
# ===========================================================================
mnn1 = np.diag(mij2000).reshape(-1, 1)
mnn0 = np.ones((TWO_NM1, 1))
mij_post = None
while ((mnn0 - mnn1) ** 2).sum() > 1e-20:
    mnn0 = mnn1.copy()
    M_mL_term = (L1base * dL[0:TWO_NM1]) / (mnn0 * L0)
    temp_local = M_mL_term.reshape(-1, 2)
    M_mL_term_ag = temp_local[:, 0:1]
    M_mL_term_na = temp_local[:, 1:2]
    Cnnjj_ag = 1 + (eta_ag + (1 - alpha) * beta_ag) * M_mL_term_ag / (alpha * beta_ag)
    Cnnjj_na = 1 + (eta_na + (1 - alpha) * beta_na) * M_mL_term_na / (alpha * beta_na)
    Cnnjj = np.column_stack([Cnnjj_ag, Cnnjj_na]).reshape(-1, 1)
    cij = Cnijk + np.eye(TWO_NM1) * np.tile(Cnnjj, (1, TWO_NM1))
    Vjmat = np.tile((dV[0:TWO_NM1] * Vi).T, (TWO_NM1, 1))
    mij_post = (cij * Vjmat) ** kappa
    mij_post = mij_post / mij_post.sum(axis=1, keepdims=True)
    mnn1 = np.mean(np.column_stack([np.diag(mij_post).reshape(-1, 1), mnn0]), axis=1, keepdims=True)


# ===========================================================================
# 16. Initial Welfare Vector (line 255)
# ===========================================================================
omega_welfare_num = (L0 / L0.sum()) * Vi * Cnnjj2000 * np.diag(mij_post).reshape(-1, 1) ** (-1.0 / kappa)
omega_welfare = omega_welfare_num / omega_welfare_num.sum()


# ===========================================================================
# 17. Base values (lines 258-261)
# ===========================================================================
mij_base = mij_post.copy()
dL_base  = dL.copy()
L1base   = L1.copy()


if __name__ == "__main__":
    print("Initial equilibrium computed.")
    print(f"  N           = {N}")
    print(f"  pi_ag shape = {pi_ag.shape}")
    print(f"  mij2000     = {mij2000.shape}")
    print(f"  Vi shape    = {Vi.shape}")
    print(f"  L0 shape    = {L0.shape}")
    print(f"  omega_welfare sum = {omega_welfare.sum():.6f}")
    print(f"  fsolve residual   = {np.max(np.abs(main_simulate(new, params))):.3e}")

# %%
