import numpy as np
# (form, year, policy, x0, t_set, wages, mean_wage, params, pm_dict_mean, bmr_dict, alert_prop, sci, pop_dict, census_grid, census_pop)
# calculate relative bad air days
def rel_ba(t, form, alert_prop, params):
    if form["baseline"] == 'none':
        r_bad = {k: v[t] for k,v in alert_prop.items()}
    return r_bad

# calculate risk
def risk(t, form, m, r_bad, params):
    R_i = {k: params["s"]*v*m[k] for k,v in r_bad.items()}
    return R_i

# calculate cost
def cost(t, form, wages, mean_wage, params, alert_prop):
    leis = params["xi"]*1*365
    if form["av_wage"] == True:
        beta_i = {k: params["EPF"]*mean_wage*leis*alert_prop[k][t] \
                                    for k,v in wages.items()}
    else:
        beta_i = {k: params["EPF"]*wages[k]*leis*alert_prop[k][t] \
                                    for k,v in wages.items()}
    return beta_i

# calculate mortalities
def morb(t, x_i, pm_dict, bmr_dict, params, eta_dict):
    dRisk = {k: params["AF"]*v[t]*0.875*0.1 for k, v in pm_dict.items()}
    return {k: params["VSL"]*params["EPF"]*(1-params["xi"]*(1-eta_dict[k][t])*v[t])*params["BMF"]*bmr_dict[k]*dRisk[k] for k, v in x_i.items()}
                                

# run full model
def model(form, year, policy, x_i, t_set, wages, mean_wage, params, pm_dict_mean, bmr_dict, alert_prop, eta_dict, sci, pop_dict, census_grid, census_pop):
    delta = {'Northeast': params["delta_n"], 'Midwest': params["delta_m"], \
            'South': params["delta_s"], 'West': params["delta_w"]}
    kappa = {'Northeast': params["kappa_n"], 'Midwest': params["kappa_m"], \
            'South': params["kappa_s"], 'West': params["kappa_w"]}

    beta = cost(t_set, form, wages, mean_wage, params, alert_prop)

    for t in t_set[:-1]:
        # calculate components of utility equation
        m_t = morb(t, x_i, pm_dict_mean, bmr_dict, params, eta_dict)
        r_bad_t = rel_ba(t, form, alert_prop, params)
        R_t = risk(t, form, m_t, r_bad_t, params)

        if t==0:
            m = {k: [v] for k, v in m_t.items()}
            r_bad = {k: [v] for k, v in r_bad_t.items()}
            R = {k: [v] for k, v in R_t.items()}
            dU = {k: [(R_t[k]-beta[k][t])/(params["EPF"]*params["R_0"])] for k in R_t.keys()}
        else:
            m = {k: m[k]+[v] for k, v in m_t.items()}
            r_bad = {k: r_bad[k]+[v] for k, v in r_bad_t.items()}
            R = {k: R[k]+[v] for k, v in R_t.items()}
            dU = {k: v+[(R_t[k]-beta[k][t])/(params["EPF"]*params["R_0"])] for k,v in dU.items()}

        delta_i = {k: sum([census_grid[k][region]*delta[region] \
                        for region in census_pop.keys()]) for k in x_i.keys()}
        kappa_i = {k: sum([census_grid[k][region]*kappa[region] \
                        for region in census_pop.keys()]) for k in x_i.keys()}

        # calculate value of social norms
        norms = {k1: sum([sci[k1][k2]*x_i[k2][t] for k2 in x_i.keys()]) for k1 in x_i.keys()}
        val_norm = {k: delta_i[k]*(2*norms[k]-1) for k in norms.keys()}

        for k, v in dU.items():
            dU[k][t] = v[t]+val_norm[k]

        if form["density"] == True:
            dx = {k: (pop_dict[k]/np.mean(list(pop_dict.values())))*kappa_i[k]\
                                *x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}
        else:
            dx = {k: kappa_i[k]*x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}

        dx = {k: kappa_i[k]*x[t]*(1-x[t])*dU[k][t] for k, x in x_i.items()}
    
        
        # update adapters for time step
        x_i = {k: x+[(max(min(x[t] + dx[k], 1), 0))] for k, x in x_i.items()}

    return x_i, dU, R, m, r_bad
