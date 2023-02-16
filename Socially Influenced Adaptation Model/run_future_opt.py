import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import os
import pickle as cp
import replicator_ms


popdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
pop = pd.read_csv(popdir + 'pop_lepeule.csv')


# Create dictionary to hold results
net_ben_dict = {}
# Set year, either 2050 or 2100
for cost_factor in [0.85]:
    net_ben_dict[cost_factor] = {}
    for year in [2050, 2100]:
        net_ben_dict[cost_factor][year] = {}
        for policy in ['REF', 'P45', 'P37']:
            net_ben_dict[cost_factor][year][policy] = {}
            for IC in ['IC' + str(i) for i in range(1,6)]:
                net_ben_dict[cost_factor][year][policy][IC] = {}
                
    
                density = 'density' # Either include population density or not ('density' or 'no_density')
    
                period = 4
                T = 5
                
                form = {'freq': 'quarterly',
                        'weights': 'exp',
                        'baseline': 'div',
                        'trend': 'annual',
                        'av_wage': True}
                
                if density == 'density':
                    form['density'] = True
                elif density == 'no_density':
                    form['density'] = False
                suffix = density
                
                indir = './input'
                outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Social Learning Model\\"
                
                # read in params from json and extract values
                with open('./params_'+suffix+'.json', 'r') as f:
                    params_raw = json.load(f)
                
                params = {key: params_raw[key]["value"] for key in params_raw \
                                                    if "value" in params_raw[key]}
                params["AF"] = (params["RR"]-1)/params["RR"]
                
                # Update parameters to make optimal
                # Social learning rate = 1
                params["kappa_m"] = 1
                params["kappa_n"] = 1
                params["kappa_s"] = 1
                params["kappa_w"] = 1
                
                params["delta_m"] = 0
                params['delta_n'] = 0
                params['delta_s'] = 0
                params['delta_w'] = 0
                
                params['b'] = params['b'] * 0
                
                params["xi"] = cost_factor
    
                
                census_pop = {'Northeast': 0.18, 'Midwest': 0.211, 'South': 0.374, 'West': 0.235}
                 
                census_df = pd.read_csv(indir+'/grid_by_census_region.csv')
                census_df.index = '(' + census_df['ROW'].astype(str) + ', ' + census_df['COL'].astype(str) + ')'
                census_df.drop(columns=['ROW', 'COL', 'geom'], inplace=True)
                census_grid = pd.DataFrame.to_dict(census_df.T)
                
                pop_df = pd.read_csv(indir+'/pop_lepeule.csv', index_col=0)
                
                wage_df = pd.read_csv(indir+'/wage_per_grid_all_b.csv')
                wage_df.index = '(' + wage_df['ROW'].astype(str) + ', ' + wage_df['COL'].astype(str) + ')'
                wages = wage_df['2000'].to_dict()
                
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
                x0 = {'Northeast': 0.2811, 'Midwest': 0.273448, \
                        'South': 0.27718, 'West': 0.474668}
                
                if year == 2050:
                    # get time series from historical fit
                    x_df = pd.read_csv(indir+'/x_fit_'+suffix+'.csv', index_col=0, parse_dates=True)
                    # get initial values from survey data to fill in missing cells
                    x0_reg = {k: [sum([census_grid[k][region]*x0[region] \
                                for region in census_pop.keys()])] for k in wages.keys()}
                    # get last time step from simulated historical data
                    x0 = x_df[-1:].to_dict(orient='list')
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
                
                pm_df = pd.read_csv(indir+'/future/'+'.'.join(['PM_quarterly', policy, IC, \
                                            'csv']), index_col=0, parse_dates=True)
                pm_df = pm_df.loc[str(start)+'-1-1':str(end)+'-12-1']
                pm_dict = pm_df.to_dict()
                pm_dict = {k: [pm for t,pm in v.items()] for k,v in pm_dict.items()}
                
                t_set = np.arange(len(pm_df))
                
                
                # run social learning model
                x, dU, R, m, BAD = replicator_ms.model(form, year, policy, x0, t_set, wages, mean_wage, params, pm_dict, bmr_dict, aqi_cur, aqi_all, sci, pop_dict, census_grid, census_pop)
                
                # convert output to dataframe
                x_df = pd.DataFrame.from_dict(x)
                x_df.index = aqi_df.index
                dU_df = pd.DataFrame.from_dict(dU)
                dU_df.index = aqi_df.index[:-1]
                R_df = pd.DataFrame.from_dict(R)
                R_df.index = aqi_df.index[:-1]
                m_df = pd.DataFrame.from_dict(m)
                m_df.index = aqi_df.index[:-1]
                bad_df = pd.DataFrame.from_dict(BAD)
                bad_df.index = aqi_df.index[:-1]
                
                # Calculate per capita cost of adaptation
                # prop adapters * hourly wage * adaptation time * proportion of time given up * days * Economic Factor
                cost_df = x_df * mean_wage * 2 * params['xi'] * 365 * params['EPF'] 
                
                # Calculate per capita all cause mortality chance
                # (Per capita economic burden of mortalities) / (VSL * EPF)
                mort_df = (m_df / (params['VSL'] * params['EPF']))
                
                out_file = '_'.join([policy, str(year), IC, suffix])
                
                # Netben for reference
                # Calculate per capita cost of adaptation
                # prop adapters * hourly wage * adaptation time * proportion of time given up * days * Economic Factor
                cost_df = x_df * mean_wage * 2 * params['xi'] * 365 * params['EPF'] 
                
                # Calculate per capita rational benefits of adaptation
                # (Economic burden under baseline adaptation) - (Economic burden under current adaptation)
                
                # ben_df = pm_df * params['AF'] * 0.1 * params['BMF'] *params['VSL'] * params['EPF'] - m_df
                
                ben_df = (m_df / (1 - params['xi']*x_df)) - m_df
                # Net Benefits = Rational benefits of adaptation - Costs of adaptation
                net_ben = ben_df - cost_df
                
                # save model output
                # x_df.to_csv(outdir+'/x_'+out_file+'.csv')
                # dU_df.to_csv(outdir+'/dU_'+out_file+'.csv')
                # R_df.to_csv(outdir+'/R_'+out_file+'.csv')
                # m_df.to_csv(outdir+'/m_'+out_file+'.csv')
                # bad_df.to_csv(outdir+'/bad_'+out_file+'.csv')
                # cost_df.to_csv(outdir+'/cost_pc_'+out_file+'.csv')
                # mort_df.to_csv(outdir+'/mort_pc_'+out_file+'.csv')
                # net_ben.to_csv(outdir+'/net_ben_'+out_file+'.csv')
                
                rc = []
                for col in net_ben:
                    rc.append(col)
                
                for cell in rc:
                    net_ben_dict[cost_factor][year][policy][IC][cell] = net_ben[cell]



rc = list(zip(pop['ROW'], pop['COL']))

tots_nb = {}
tots_b = {}
tots_c = {}
for cost_factor in [0.85]:
    tots_nb[cost_factor] = {}
    tots_b[cost_factor] = {}
    tots_c[cost_factor] = {}
    for year in [2050, 2100]:
        tots_nb[cost_factor][year] = {}
        tots_b[cost_factor][year] = {}
        tots_c[cost_factor][year] = {}
        for policy in ['REF', 'P45', 'P37']:
            tots_nb[cost_factor][year][policy] = []
            tots_b[cost_factor][year][policy] = []
            tots_c[cost_factor][year][policy] = []
            for IC in ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']:
                totlist_nb = []
                totlist_b = []
                totlist_c = []
                for r, c in rc:
                    popu = pop.loc[(pop['ROW'] == r)&(pop['COL'] == c), str(year)].item()
                    totlist_nb.append(np.mean(net_ben_dict[cost_factor][year][policy][IC][str((r, c))]) * popu)
                    
                
                tots_nb[cost_factor][year][policy].append(sum(totlist_nb))
                
            
            tots_nb[cost_factor][year][policy] = sum(tots_nb[cost_factor][year][policy]) / len(tots_nb[cost_factor][year][policy]) / 1E12
            
            







            
# with open(outdir + 'net_ben_opt', 'wb') as f:
#     cp.dump(net_ben_dict,f)
#     f.close()           
                       