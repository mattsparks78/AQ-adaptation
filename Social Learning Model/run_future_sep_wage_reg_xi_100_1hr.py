import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import sys
import os
import pickle as cp
import replicator_av_wage_reg_xi_100_1hr
import copy

inf_rate = 0.2

# Set conversion factor for 2000 to 2005 USD for wage
w_conv = 1.13


# Create dict to hold ben, cost, net_ben values for each case
ben_dict = {}
cost_dict = {}
net_ben_dict = {}
# Set year, either 2050 or 2100
for k_factor in [1]:
    ben_dict[k_factor] = {}
    cost_dict[k_factor] = {}
    net_ben_dict[k_factor] = {}
    for cost_factor in [1]:
        ben_dict[k_factor][cost_factor] = {}
        cost_dict[k_factor][cost_factor] = {}
        net_ben_dict[k_factor][cost_factor] = {}
        for risk_factor in [1]:
            ben_dict[k_factor][cost_factor][risk_factor] = {}
            cost_dict[k_factor][cost_factor][risk_factor] = {}
            net_ben_dict[k_factor][cost_factor][risk_factor] = {}
            for year in [2050, 2100]:
                ben_dict[k_factor][cost_factor][risk_factor][year] = {}
                cost_dict[k_factor][cost_factor][risk_factor][year] = {}
                net_ben_dict[k_factor][cost_factor][risk_factor][year] = {}
                for policy in ['REF', 'P45', 'P37']:
                    ben_dict[k_factor][cost_factor][risk_factor][year][policy] = {}
                    cost_dict[k_factor][cost_factor][risk_factor][year][policy] = {}
                    net_ben_dict[k_factor][cost_factor][risk_factor][year][policy] = {}
                    for IC in ['IC' + str(i) for i in range(1,6)]:
                        ben_dict[k_factor][cost_factor][risk_factor][year][policy][IC] = {}
                        cost_dict[k_factor][cost_factor][risk_factor][year][policy][IC] = {}
                        net_ben_dict[k_factor][cost_factor][risk_factor][year][policy][IC] = {}
                        density = 'density' # Either include population density or not ('density' or 'no_density')
            
                        period = 4
                        T = 5
                        
                        form = {'freq': 'quarterly',
                                'weights': 'exp',
                                'baseline': 'none',
                                'trend': 'annual',
                                'av_wage': False}
                        
                        if density == 'density':
                            form['density'] = True
                        elif density == 'no_density':
                            form['density'] = False
                        suffix = 'density_sep_wage_regional_xi'
                        
                        indir = './input'
                        indir1 = r'.\Data\Future AQ\\'
                        outdir = r".\Data\Social Learning Model\Sep 100\\"
                        
                        # read in params from json and extract values
                        with open('./params_'+suffix+'.json', 'r') as f:
                            params_raw = json.load(f)
                        
                        params = {key: params_raw[key]["value"] for key in params_raw \
                                                            if "value" in params_raw[key]}
                        params["AF"] = (params["RR"]-1)/params["RR"]
                        
                        # Change parameters that need changing
                        params['kappa_m'] = params['kappa_m'] * k_factor
                        params['kappa_n'] = params['kappa_n'] * k_factor
                        params['kappa_s'] = params['kappa_s'] * k_factor
                        params['kappa_w'] = params['kappa_w'] * k_factor
                        
                        
                        census_pop = {'Northeast': 0.18, 'Midwest': 0.211, 'South': 0.374, 'West': 0.235}
                         
                        census_df = pd.read_csv(indir+'/grid_by_census_region.csv')
                        census_df.index = '(' + census_df['ROW'].astype(str) + ', ' + census_df['COL'].astype(str) + ')'
                        census_df.drop(columns=['ROW', 'COL', 'geom'], inplace=True)
                        census_grid = pd.DataFrame.to_dict(census_df.T)
                        
                        pop_df = pd.read_csv(indir+'/pop_lepeule.csv', index_col=0)
                        
                        wage_df = pd.read_csv(indir+'/wage_per_grid_all_b.csv')
                        wage_df['2005'] = wage_df['2000'] * w_conv * cost_factor
                        wage_df.index = '(' + wage_df['ROW'].astype(str) + ', ' + wage_df['COL'].astype(str) + ')'
                        wages = wage_df['2005'].to_dict()
                        
                        # delete grid cells not found in other input files
                        del wages['(1, 20)']
                        del wages['(14, 24)']
                        del wages['(14, 25)']
                        wage_list = list(wages.values())
                        
                        sci = pd.read_pickle(indir+'/sci_normalized_b')
                        sci = {str(k1): {str(k2): v2 for k2, v2 in sci[k1].items()} for k1, v1 in sci.items()}
                        
                        bmr = pd.read_csv(indir+'/bmr.csv', index_col=0, parse_dates=True)
                        bmr_dict = bmr['BaseMortRate'].to_dict()
                        
                        # census region x0's from PN survey - used for grid cells with missing data in histirical runs
                        x0_region = {'Northeast': 0.2811, 'Midwest': 0.273448, \
                                'South': 0.27718, 'West': 0.474668}
                        
                        if year == 2050:
                            # get time series from historical fit
                            x_df = pd.read_csv(indir+'/x_0_'+suffix+'.csv')
                            # get initial values from survey data to fill in missing cells
                            x0_reg = {k: [sum([census_grid[k][region]*x0_region[region] \
                                        for region in census_pop.keys()])] for k in wages.keys()}
                            # get first time step from simulated historical data
                            x0 = x_df.iloc[0].to_dict()
                            for k in x0.keys():
                                x0[k] = [x0[k]]
                            for k, v in x0_reg.items():
                                # fill in missing cells with PN survey data
                                if k not in x0: x0[k] = v
                        elif year == 2100:
                            # get time series from 2050 run
                            in_file = '_'.join([policy, str(year-50), IC, suffix])
                            x_df = pd.read_csv(outdir+'/x_'+in_file+'.csv', index_col=0, parse_dates=True)
                            x0 = x_df[-1:].to_dict(orient='list')
                        
                        # update params for year and policy
                        params["EPF"] = params_raw["EPF"][str(year)][policy]
                        params["BMF"] = params_raw["BMF"][str(year)]
                        start = year-9
                        end = year+15
                        
                        pop_dict = pop_df[str(year)].to_dict()
                        mean_wage = sum([wages[k]*(pop_dict[k]/sum(pop_dict.values())) for k in pop_dict.keys()]) 
                        reg_dict = {region: {k: census_grid[k][region]*pop_dict[k] \
                                            for k in pop_dict.keys()} for region in census_pop.keys()}
                        
                        aqi_df = pd.read_csv(indir+'/future/'+'.'.join(['BAD', policy, IC, 'csv']), \
                                                                index_col=0, parse_dates=True)
                        aqi_df.drop(columns=['prop_total'], inplace=True)
                        
                        if form["trend"] == 'annual':
                            aqi_df = aqi_df.rolling(period).mean()
                        aqi_df = aqi_df.fillna(method="bfill")
                        
                        aqi_df = aqi_df.loc[str(start-5)+'-1-1':str(end)+'-12-1']
                        aqi_all = {}
                        aqi_cur = {}
                        
                        for column in aqi_df:
                            aqi_all[column] = np.array(aqi_df[column])
                            aqi_cur[column] = aqi_all[column][T*period:]
                            aqi_diff = len(aqi_all[column])-len(aqi_cur[column])
                        
                        aqi_df = aqi_df.loc[str(start)+'-1-1':str(end)+'-12-1']
                        
                        fname = policy + '_PM25_2050_2100_d'
                        with open(indir1 + fname, 'rb') as f:
                             test = cp.load(f)
                        f.close()
                        
                        pm_dict = {k: [] for k in pop_dict}
                        if year == 2050:
                            year_set = range(2041,2066)
                        elif year == 2100:
                            year_set = range(2091, 2116)
                        
                        del test['2041']['IC1'][(14,24)]
                        
                        
                        for cell in test['2041']['IC1'].keys():
                            for Yr in year_set:
                                for q in ['q1','q2','q3','q4']:
                                    if q == 'q1':
                                        pm_dict[str(cell)].append(np.concatenate((test[str(Yr-1)][IC][cell][range(90,365)], test[str(Yr)][IC][cell][range(0,90)]), axis = 0))
                                    elif q == 'q2':
                                        pm_dict[str(cell)].append(np.concatenate((test[str(Yr-1)][IC][cell][range(181,365)], test[str(Yr)][IC][cell][range(0,181)]), axis = 0))
                                    elif q == 'q3':
                                        pm_dict[str(cell)].append(np.concatenate((test[str(Yr-1)][IC][cell][range(273,365)], test[str(Yr)][IC][cell][range(0,273)]), axis = 0))
                                    elif q == 'q4':
                                        pm_dict[str(cell)].append(test[str(Yr)][IC][cell][range(0,365)])
                        
                        
                        pm_df = pd.read_csv(indir+'/future/'+'.'.join(['PM_quarterly', policy, IC, \
                                                    'csv']), index_col=0, parse_dates=True)
                        pm_df = pm_df.loc[str(start)+'-1-1':str(end)+'-12-1']
                        pm_df = pm_df.drop(['(14, 24)'], axis = 1)
                        pm_dict_mean = pm_df.to_dict()
                        pm_dict_mean = {k: [pm for t,pm in v.items()] for k,v in pm_dict_mean.items()}
                        
                        t_set = np.arange(len(pm_df))
                        
                        alert_prop = {k: [] for k in pop_dict}
                        for cell in alert_prop.keys():
                            for i in t_set:
                                year_pm = pm_dict[str(cell)][i]
                                num_bad = len(year_pm[year_pm > 35.5])
                                alert_prop[cell].append(num_bad / 365)
                            alert_prop[cell] = np.array(alert_prop[cell])
                        
                        # Calculate psi by region
                        t_out_base = 1
                        t_in_base = 24 - 1
                        
                        psi_reg = {}
                        for reg in ['Northeast', 'South', 'Midwest', 'West']:
                            if reg == 'Northeast':
                                xi_temp = params['xi_n']
                            elif reg == 'South':
                                xi_temp = params['xi_s']
                            elif reg == 'Midwest':
                                xi_temp = params['xi_m']
                            elif reg == 'West':
                                xi_temp = params['xi_w']    
                            t_out_ad = (1 - xi_temp) * t_out_base
                            t_in_ad = 24 - t_out_ad
                            exp_ad = t_in_ad * inf_rate + t_out_ad
                            exp_base = t_in_base * inf_rate + t_out_base
                            psi_reg[reg] = exp_ad / exp_base
                    
                        xi_reg = {'Northeast': params['xi_n'], 'South':params['xi_s'], \
                                  'Midwest':params['xi_m'], 'West': params['xi_w']}
                        
                        # Create xi dict
                        xi_dict = {k: 0 for k in pop_dict}
                        for k in pop_dict.keys():
                            for reg in ['Northeast', 'Midwest', 'South', 'West']:
                                xi_dict[k] = xi_dict[k] + xi_reg[reg] * census_grid[k][reg]
                        
                    
                        # Calculate psi by cell based on regional proportion
                        psi_cell = {k: 0 for k in pop_dict}
                        for k in psi_cell:
                            for reg in ['Northeast', 'Midwest', 'South', 'West']:
                                psi_cell[k] = psi_cell[k] + psi_reg[reg] * census_grid[k][reg]
                        
                        #Calculate eta
                        eta_dict = {k: [] for k in pop_dict}
                        for cell in alert_prop.keys():
                            for i in t_set:
                                year_pm = pm_dict[str(cell)][i]
                                alert_days = year_pm[year_pm >= 35.5]
                                nonalert_days = year_pm[year_pm < 35.5]
                                alert_total = sum(alert_days) * psi_cell[cell]
                                nonalert_total = sum(nonalert_days)
                                eta_dict[cell].append((alert_total + nonalert_total) / sum(year_pm))
                   
                        # run social learning model
                        x, dU, R, m, BAD = replicator_av_wage_reg_xi_100_1hr.model(form, year, policy, x0, t_set, wages, mean_wage, params, pm_dict_mean, bmr_dict, alert_prop, eta_dict, sci, pop_dict, census_grid, census_pop)
                        
                        # convert output to dataframe
                        x_df = pd.DataFrame.from_dict(x)
                        x_df = x_df.drop(99, axis = 0)
                        x_df.index = aqi_df.index[:-1]
                        dU_df = pd.DataFrame.from_dict(dU)
                        dU_df.index = aqi_df.index[:-1]
                        R_df = pd.DataFrame.from_dict(R)
                        R_df.index = aqi_df.index[:-1]
                        m_df = pd.DataFrame.from_dict(m)
                        m_df.index = aqi_df.index[:-1]
                        bad_df = pd.DataFrame.from_dict(BAD)
                        bad_df.index = aqi_df.index[:-1]
                        
                        
                        
                        # Calculate per capita all cause mortality chance
                        # (Per capita economic burden of mortalities) / (VSL * EPF)
                        mort_df = (m_df / (params['VSL'] * params['EPF']))
                        
                        out_file = '_'.join([policy, str(year), IC, suffix])
                        
                        # Netben for reference
                        # Calculate per capita cost of adaptation
                        # prop adapters * hourly wage * adaptation time * proportion of time given up * days * Economic Factor
                        alert_df = pd.DataFrame.from_dict(alert_prop)
                        alert_df = alert_df.drop(99, axis = 0)
                        alert_df.index = aqi_df.index[:-1]
                        
                        
                        # Calculate per capita cost of 
                        cost_df = copy.deepcopy(x_df)
                        
                        for cell in alert_prop.keys():
                            cost_df[cell] = cost_df[cell] * wages[cell] * xi_dict[cell] * 365 * params['EPF'] * alert_df[cell]
                            
                                                
                        # Calculate per capita rational benefits of adaptation
                        # (Economic burden under baseline adaptation) - (Economic burden under current adaptation)
                        eta_df = pd.DataFrame.from_dict(eta_dict)
                        eta_df = eta_df.drop(99, axis = 0)
                        eta_df.index = aqi_df.index[:-1]
                        ben_df = (m_df / (1 - x_df*(1-eta_df))) - m_df
                        # Net Benefits = Rational benefits of adaptation - Costs of adaptation
                        net_ben = ben_df - cost_df
                        
                        rc = []
                        
                        for col in net_ben:
                            rc.append(col)
                        
                        for cell in rc:
                            net_ben_dict[k_factor][cost_factor][risk_factor][year][policy][IC][cell] = net_ben[cell]
                            ben_dict[k_factor][cost_factor][risk_factor][year][policy][IC][cell] = ben_df[cell]
                            cost_dict[k_factor][cost_factor][risk_factor][year][policy][IC][cell] = cost_df[cell]
                        
                        
                        # save model output
                        x_df.to_csv(outdir+'/x_'+out_file+'.csv')
                        dU_df.to_csv(outdir+'/dU_'+out_file+'.csv')
                        R_df.to_csv(outdir+'/R_'+out_file+'.csv')
                        m_df.to_csv(outdir+'/m_'+out_file+'.csv')
                        bad_df.to_csv(outdir+'/bad_'+out_file+'.csv')
                        cost_df.to_csv(outdir+'/cost_pc_'+out_file+'.csv')
                        mort_df.to_csv(outdir+'/mort_pc_'+out_file+'.csv')
                        net_ben.to_csv(outdir+'/net_ben_'+out_file+'.csv')

# with open(outdir + 'net_ben_dict_sep_wage_reg_xi_100_1hr', 'wb') as f:
#     cp.dump(net_ben_dict,f)
#     f.close()                           
            
# with open(outdir + 'cost_dict_sep_wage_reg_xi_100_1hr', 'wb') as f:
#     cp.dump(cost_dict,f)
#     f.close()    

# with open(outdir + 'ben_dict_sep_wage_reg_xi_100_1hr', 'wb') as f:
#     cp.dump(ben_dict,f)
#     f.close()                 