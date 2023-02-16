import numpy as np

# calculate moving baseline of bad air days


def mvg_bl(t, form, aqi_diff, aqi_all, params):
    B = 0
    if form["weights"] == 'lin':
        w = [max(params["b"]-(params["a"]*k), 0) for k in range(params["T"])]
    elif form["weights"] == 'exp':
        w = [np.exp(-params["a"]*k) for k in range(params["T"])]
    w_norm = [w[k]/sum(w) for k in range(params["T"])]
    t_del = np.array([t+(aqi_diff-1)-(params["period"]*k)
                     for k in range(params["T"])])
    aqi_t_del = {k: v[t_del] for k, v in aqi_all.items()}
    B = {k: sum([v[i]*w_norm[i] for i in range(params["T"])])
         for k, v in aqi_t_del.items()}
    return B

# calculate relative bad air days


def rel_ba(t, form, aqi_cur, aqi_all, params):
    arb_grid = list(aqi_all.keys())[0]
    aqi_diff = len(aqi_all[arb_grid])-len(aqi_cur[arb_grid])
    bl = mvg_bl(t, form, aqi_diff, aqi_all, params)
    if form["baseline"] == 'sub':
        r_bad = {k: v[t]-bl[k] for k, v in aqi_cur.items()}
    elif form["baseline"] == 'div':
        r_bad = {k: (v[t]+params["eps"])/(bl[k]+params["eps"])
                 for k, v in aqi_cur.items()}
    return r_bad

# calculate risk


def risk(t, form, m, r_bad, aqi_cur, aqi_all, params):
    R_i = {k: params["s"]*v*m[k]+params["EPF"]*params["R_0"]
           for k, v in r_bad.items()}
    return R_i

# calculate cost


def cost(t, form, wages, mean_wage, params):
    leis = params["xi"]*2*365
    if form["av_wage"] == True:
        beta_i = {k: params["EPF"]*mean_wage*leis
                  for k, v in wages.items()}
    else:
        beta_i = {k: params["EPF"]*wages[k]*leis
                  for k, v in wages.items()}
    return beta_i

# calculate mortalities


def morb(t, x_i, pm_dict, bmr_dict, params):
    """
    

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    x_i : TYPE
        DESCRIPTION.
    pm_dict : TYPE
        DESCRIPTION.
    bmr_dict : TYPE
        DESCRIPTION.
    params : TYPE
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    dRisk = {k: params["AF"]*v[t]*0.1 for k, v in pm_dict.items()}
    return {k: params["VSL"]*params["EPF"]*(1-params["xi"]*v[t])*params["BMF"]
            * bmr_dict[k]*dRisk[k] for k, v in x_i.items()}

# run full model


def model(form, year, policy, x_i, t_set, wages, mean_wage, params, pm_dict, bmr_dict, aqi_cur, aqi_all, sci, pop_dict, census_grid, census_pop):

    delta = {'Northeast': params["delta_n"], 'Midwest': params["delta_m"],
             'South': params["delta_s"], 'West': params["delta_w"]}
    kappa = {'Northeast': params["kappa_n"], 'Midwest': params["kappa_m"],
             'South': params["kappa_s"], 'West': params["kappa_w"]}

    beta = cost(t_set, form, wages, mean_wage, params)

    for t in t_set[:-1]:
        # calculate components of utility equation
        m_t = morb(t, x_i, pm_dict, bmr_dict, params)
        r_bad_t = rel_ba(t, form, aqi_cur, aqi_all, params)
        R_t = risk(t, form, m_t, r_bad_t, aqi_cur, aqi_all, params)

        if t == 0:
            m = {k: [v] for k, v in m_t.items()}
            r_bad = {k: [v] for k, v in r_bad_t.items()}
            R = {k: [v] for k, v in R_t.items()}
            dU = {k: [(R_t[k]-beta[k])/(params["EPF"]*params["R_0"])]
                  for k in R_t.keys()}
        else:
            m = {k: m[k]+[v] for k, v in m_t.items()}
            r_bad = {k: r_bad[k]+[v] for k, v in r_bad_t.items()}
            R = {k: R[k]+[v] for k, v in R_t.items()}
            dU = {k: v+[(R_t[k]-beta[k])/(params["EPF"]*params["R_0"])]
                  for k, v in dU.items()}

        delta_i = {k: sum([census_grid[k][region]*delta[region]
                           for region in census_pop.keys()]) for k in x_i.keys()}
        kappa_i = {k: sum([census_grid[k][region]*kappa[region]
                           for region in census_pop.keys()]) for k in x_i.keys()}

        # calculate value of social norms
        norms = {k1: sum([sci[k1][k2]*x_i[k2][t]
                         for k2 in x_i.keys()]) for k1 in x_i.keys()}
        val_norm = {k: delta_i[k]*(2*norms[k]-1) for k in norms.keys()}

        for k, v in dU.items():
            dU[k][t] = v[t]+val_norm[k]

        if form["density"] == True:
            dx = {k: (pop_dict[k]/np.mean(list(pop_dict.values())))*kappa_i[k]
                  * x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}
        else:
            dx = {k: kappa_i[k]*x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}

        dx = {k: kappa_i[k]*x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}
        # update adapters for time step
        x_i = {k: x+[(max(min(x[t] + dx[k], 1), 0))] for k, x in x_i.items()}

    return x_i, dU, R, m, r_bad
