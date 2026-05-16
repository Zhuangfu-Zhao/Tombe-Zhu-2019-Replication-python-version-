import numpy as np

def eqm(X, params):

    N       = params["N"]           # number of regions 
    L0      = params["L0"]          # registration labor distribution
    La      = params["La"]          # labor distribution in agriculture after migration
    Ln      = params["Ln"]          # labor distribution in non-agriculture after migration
    L1base  = params["L1base"]      # labor distribution

    Ir       = params["Ir"]         # total income (rural)
    Iu       = params["Iu"]         # total income (urban)
    R_ag     = params["R_ag"]       # revenue of agriculture
    R_na     = params["R_na"]       # revenue of non-agriculture

    # production params
    theta   = params["theta"]       # dispersion of distribution of productivity
    beta_ag = params["beta_ag"]     # labor share of agriculture production
    beta_na = params["beta_na"]     # labor share of non-agriculture production
    eta_ag  = params["eta_ag"]      # land share of agriculture production
    eta_na  = params["eta_na"]      # labor share of non-agriculture production
    sigma   = params["sigma"]       # intermediate share matrix

    # preference params
    kappa   = params["kappa"]       # dispersion of distribution of preference
    alpha   = params["alpha"]       # goods share (ag + na) of preference
    psi     = params["psi"]         # ag share in goods

    # trade
    pi_ag   = params["pi_ag"]       # trade shares (agriculture)
    pi_na   = params["pi_na"]       # trade shares (non-agriculture)

    # migration
    Vi      = params["Vi"]          # real income per worker
    mij2000 = params["mij2000"]     # migration share relative to registration 
    Cnijk   = params["Cnijk"]       # migration cost

    # deviations 
    dTa      = params["dTa"]
    dTn      = params["dTn"]
    dni_ag   = params["dni_ag"]
    dni_na   = params["dni_na"] 


    # eqm components
    X = X.reshape(-1, 1)
    dwa = X[0:N]
    dwn = X[N:2 * N]
    dPa = X[2 * N:3 * N]
    dPn = X[3 * N:4 * N]
    dLa = X[4 * N:5 * N]
    dLn = X[5 * N:6 * N]
    dL = np.column_stack([dLa, dLn]).reshape(-1, 1)

    # changes in land rents
    dra = dwa * dLa 
    drn = dwn * dLn

    # new sectoral revenues
    R_ag_new = dwa * dLa * R_ag
    R_na_new = dwn * dLn * R_na

    # new household final demand
    Ir_new = (beta_ag + eta_ag) * R_ag_new / alpha
    dIr = Ir_new / Ir
    Iu_new = (beta_na + eta_na) * R_na_new / alpha
    dIu = Iu_new / Iu
    Da_new_r = alpha * psi * Ir_new
    Da_new_u = alpha * psi * Iu_new
    Dn_new_r = (alpha * (1 - psi)) * Ir_new
    Dn_new_u = (alpha * (1 - psi)) * Iu_new

    # new trade shares
    inside_ag = (dwa.T ** beta_ag
                 * (dPa.T ** sigma[0, 0] * dPn.T ** sigma[0, 1]) ** (1 - beta_ag - eta_ag)
                 * dra.T ** eta_ag
                 / (dTa.T ** (1.0 / theta)))
    kron_inside_ag = np.kron(inside_ag, np.ones((N, 1)))
    denominator_ag = np.sum(pi_ag * (dni_ag * kron_inside_ag) ** (-theta), axis=1, keepdims=True)
    pi_ag_new = (pi_ag * (dni_ag * kron_inside_ag) ** (-theta)) / np.kron(denominator_ag, np.ones((1, N)))

    inside_na = (dwn.T ** beta_na
                 * (dPa.T ** sigma[1, 0] * dPn.T ** sigma[1, 1]) ** (1 - beta_na - eta_na)
                 * drn.T ** eta_na
                 / (dTn.T ** (1.0 / theta)))
    kron_inside_na = np.kron(inside_na, np.ones((N, 1)))
    denominator_na = np.sum(pi_na * (dni_na * kron_inside_na) ** (-theta), axis=1, keepdims=True)
    pi_na_new = (pi_na * (dni_na * kron_inside_na) ** (-theta)) / np.kron(denominator_na, np.ones((1, N)))

    # new sector expenditures
    X_ag_new = np.linalg.solve(pi_ag_new.T, R_ag_new)
    X_na_new = np.linalg.solve(pi_na_new.T, R_na_new)

    # real income per worker
    dP = dPa ** psi * dPn ** (1 - psi)
    dVa = dwa / (dP ** alpha * dra ** (1 - alpha))
    dVn = dwn / (dP ** alpha * drn ** (1 - alpha))
    dV = np.column_stack([dVa, dVn]).reshape(-1, 1)

    # migration shares (inner fixed-point loop)
    NN = 2 * (N - 1)
    mnn1 = np.diag(mij2000).reshape(-1, 1)
    mnn0 = np.ones((NN, 1))
    mij_loc = None

    while ((mnn0 - mnn1) ** 2).sum() > 1e-20:
        mnn0 = mnn1.copy()
        M_mL_term = (L1base * dL[0:NN]) / (mnn0 * L0)
        temp_local = M_mL_term.reshape(-1, 2)
        M_mL_term_ag_l = temp_local[:, 0:1]
        M_mL_term_na_l = temp_local[:, 1:2]
        Cnnjj_ag_l = 1 + (eta_ag + (1 - alpha) * beta_ag) * M_mL_term_ag_l / (alpha * beta_ag)
        Cnnjj_na_l = 1 + (eta_na + (1 - alpha) * beta_na) * M_mL_term_na_l / (alpha * beta_na)
        Cnnjj_l = np.column_stack([Cnnjj_ag_l, Cnnjj_na_l]).reshape(-1, 1)
        cij_l = Cnijk + np.eye(NN) * np.tile(Cnnjj_l, (1, NN))
        Vjmat_l = np.tile((dV[0:NN] * Vi).T, (NN, 1))
        mij_loc = (cij_l * Vjmat_l) ** kappa
        mij_loc = mij_loc / mij_loc.sum(axis=1, keepdims=True)
        mnn1 = np.mean(np.column_stack([np.diag(mij_loc).reshape(-1, 1), mnn0]), axis=1, keepdims=True)

    # new labour allocation
    L1_loc = (L0.T @ mij_loc).T

    # system of equations
    F = np.vstack([
        dPa - denominator_ag ** (-1.0 / theta),
        dPn - denominator_na ** (-1.0 / theta),
        L1_loc[0:NN] / L1base[0:NN] - dL[0:NN],
        np.array([[(dLa[N - 1, 0] * La[N - 1, 0] + dLn[N - 1, 0] * Ln[N - 1, 0]
                    - (La[N - 1, 0] + Ln[N - 1, 0]))]]),
        np.array([[dwa[N - 1, 0] - dwn[N - 1, 0]]]),
        X_ag_new - (Da_new_r + Da_new_u
                    + (1 - beta_ag - eta_ag) * sigma[0, 0] * R_ag_new
                    + (1 - beta_na - eta_na) * sigma[1, 0] * R_na_new),
        X_na_new - (Dn_new_r + Dn_new_u
                    + (1 - beta_na - eta_na) * sigma[1, 1] * R_na_new
                    + (1 - beta_ag - eta_ag) * sigma[0, 1] * R_ag_new),
    ])
    return F.flatten()