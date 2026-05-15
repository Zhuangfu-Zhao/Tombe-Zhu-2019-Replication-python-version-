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
    sigma_p    = params["sigma"]
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
                 * (dPa.T ** sigma_p[0, 0] * dPn.T ** sigma_p[0, 1]) ** (1 - beta_ag_p - eta_ag_p)
                 * dra.T ** eta_ag_p
                 / (dTa_p.T ** (1.0 / theta_p)))
    kron_inside_ag = np.kron(inside_ag, np.ones((Np, 1)))
    denominator_ag = np.sum(pi_ag_p * (dni_ag_p * kron_inside_ag) ** (-theta_p), axis=1, keepdims=True)
    pi_ag_new = (pi_ag_p * (dni_ag_p * kron_inside_ag) ** (-theta_p)) / np.kron(denominator_ag, np.ones((1, Np)))

    inside_na = (dwn.T ** beta_na_p
                 * (dPa.T ** sigma_p[1, 0] * dPn.T ** sigma_p[1, 1]) ** (1 - beta_na_p - eta_na_p)
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
    NN = 2 * (Np - 1)
    mnn1 = np.diag(mij2000_p).reshape(-1, 1)
    mnn0 = np.ones((NN, 1))
    mij_loc = None
    while ((mnn0 - mnn1) ** 2).sum() > 1e-20:
        mnn0 = mnn1.copy()
        M_mL_term = (L1base_p * dL[0:NN]) / (mnn0 * L0_p)
        temp_local = M_mL_term.reshape(-1, 2)
        M_mL_term_ag_l = temp_local[:, 0:1]
        M_mL_term_na_l = temp_local[:, 1:2]
        Cnnjj_ag_l = 1 + (eta_ag_p + (1 - alpha_p) * beta_ag_p) * M_mL_term_ag_l / (alpha_p * beta_ag_p)
        Cnnjj_na_l = 1 + (eta_na_p + (1 - alpha_p) * beta_na_p) * M_mL_term_na_l / (alpha_p * beta_na_p)
        Cnnjj_l = np.column_stack([Cnnjj_ag_l, Cnnjj_na_l]).reshape(-1, 1)
        cij_l = Cnijk_p + np.eye(NN) * np.tile(Cnnjj_l, (1, NN))
        Vjmat_l = np.tile((dV[0:NN] * Vi_p).T, (NN, 1))
        mij_loc = (cij_l * Vjmat_l) ** kappa_p
        mij_loc = mij_loc / mij_loc.sum(axis=1, keepdims=True)
        mnn1 = np.mean(np.column_stack([np.diag(mij_loc).reshape(-1, 1), mnn0]), axis=1, keepdims=True)

    # New labour allocation
    L1_loc = (L0_p.T @ mij_loc).T

    # System of equations
    F = np.vstack([
        dPa - denominator_ag ** (-1.0 / theta_p),
        dPn - denominator_na ** (-1.0 / theta_p),
        L1_loc[0:NN] / L1base_p[0:NN] - dL[0:NN],
        np.array([[(dLa[Np - 1, 0] * La_data_p[Np - 1, 0] + dLn[Np - 1, 0] * Ln_data_p[Np - 1, 0]
                    - (La_data_p[Np - 1, 0] + Ln_data_p[Np - 1, 0]))]]),
        np.array([[dwa[Np - 1, 0] - dwn[Np - 1, 0]]]),
        X_ag_new - (Da_new_r + Da_new_u
                    + (1 - beta_ag_p - eta_ag_p) * sigma_p[0, 0] * R_ag_new
                    + (1 - beta_na_p - eta_na_p) * sigma_p[1, 0] * R_na_new),
        X_na_new - (Dn_new_r + Dn_new_u
                    + (1 - beta_na_p - eta_na_p) * sigma_p[1, 1] * R_na_new
                    + (1 - beta_ag_p - eta_ag_p) * sigma_p[0, 1] * R_ag_new),
    ])
    return F.flatten()