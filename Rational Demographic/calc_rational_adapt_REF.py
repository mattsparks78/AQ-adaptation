# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:09:50 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
from shapely import wkt

# Set directories
indir1 = r".\Data\General\\"
indir2 = r".\Data\Future AQ\\"
indir3 = r".\Data\AGU2020-Fig3\\"
outdir = r".\Data\Rational Actor Infiltration Dem\\"

# Load future concentration data. REF, P3.7 or P4.5
with open(indir2 + 'REF_PM25_2050_2100_d', 'rb') as f:
    pm_out = cp.load(f)
f.close()

# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_race.csv")

# Import population data
pop = pd.read_csv(indir1 + 'pop_lepeule_race.csv')

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']), list(wage['COL'])))

# Set relative risk 
RR = 1.14

# Calculate attributable fraction
AF = (RR - 1) / RR

# Set VSL. Converted from $2006 USD to 2005
VSL = 7104000

# Set conversion factor for 2000 to 2005 USD for wage
w_conv = 1.13

# Proportion of population who adapt. Assumed 1 here. Can update with replicator numbers.
ppa = 1
# Calculate fractions of indoor and outdoor time
t_out = 2
t_in = 24 - t_out

prop_out = t_out / 24
prop_in = t_in / 24

# Set adaptation time in hours. 
AT = 2

# Set infiltration rate based on empirical data/for sensitvity analysis
inf_rate = 0.2

# Calculate ratio of adaptation exposure to baseline exposure
exp_base = prop_out + prop_in*inf_rate
exp_ad = (prop_out + prop_in)*inf_rate

ad_reduct = exp_ad / exp_base
#%%

# Calculate adaptation days for all scenarios
# Initialize dictionaries
pm_out_base = {}
pm_out_adapt = {}
adapt_count = {}
adapt_cost = {}
adapt_ben = {}
adapt_util = {}
adapt_cost_pc = {}
adapt_ben_pc = {}
adapt_util_pc = {}

# Iterate over both sets of years (for 2050 and 2100)
for years in year_set:
    # Iterate over all years in year set
    if np.mean(years) < 2070:
        # Economic projection factor
        EPF = 2.01
        # Population projection factor
        PPF = 1.48
        # PM Mortality Factor (for baseline mortaltiy)
        PMMF = 1.12
        # Set VSL. Converted from $2006 USD to 2005
        VSL = 7104000
        # Update VSL based on EPF
        VSL = VSL * EPF
    else:
        EPF = 3.94
        PPF = 1.73
        PMMF = 0.92
        # Set VSL. Converted from $2006 USD to 2005
        VSL = 7104000
        # Update VSL based on EPF
        VSL = VSL * EPF
    for Yr in years:
        # Initialize next level of dictionaries
        pm_out_base[Yr] = {}
        pm_out_adapt[Yr] = {}
        adapt_count[Yr] = {}
        adapt_cost[Yr] = {}
        adapt_ben[Yr] = {}
        adapt_util[Yr] = {}
        adapt_cost_pc[Yr] = {}
        adapt_ben_pc[Yr] = {}
        adapt_util_pc[Yr] = {}

        # Iterate over all initial conditions in each year
        for IC in IC_list:
            # Initialize next level of dictionaries
            pm_out_base[Yr][IC] = {}
            pm_out_adapt[Yr][IC] = {}
            adapt_count[Yr][IC] = {}
            adapt_cost[Yr][IC] = {}
            adapt_ben[Yr][IC] = {}
            adapt_util[Yr][IC] = {}
            adapt_cost_pc[Yr][IC] = {}
            adapt_ben_pc[Yr][IC] = {}
            adapt_util_pc[Yr][IC] = {}

            # Iterate over all cells in each scenario
            for r, c in rc:
                # Initialize next level of dictionaries
                pm_out_base[Yr][IC][(r, c)] = {}
                pm_out_adapt[Yr][IC][(r, c)] = {}
                adapt_count[Yr][IC][(r, c)] = {}
                adapt_cost[Yr][IC][(r, c)] = {}
                adapt_ben[Yr][IC][(r, c)] = {}
                adapt_util[Yr][IC][(r, c)] = {}
                adapt_cost_pc[Yr][IC][(r, c)] = {}
                adapt_ben_pc[Yr][IC][(r, c)] = {}
                adapt_util_pc[Yr][IC][(r, c)] = {}
                
                # Calculate baseline mean outdoor PM concentration
                pm_o = np.mean(pm_out[str(Yr)][IC][(r, c)])
                # Store baseline mean outdoor PM concentration in dictionary
                pm_out_base[Yr][IC][(r, c)] = pm_o
                # Find and store baseline mortality rate for the current grid cell
                base_mort = bm.loc[(bm['ROW'] == r) & (
                    bm['COL'] == c), 'BaseMortRate'].item() * PMMF
                
                # Iterate over each racial demographic
                for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
               
                    
                    # Find and store mean hourly wage for the current grid cell
                    w_i = wage.loc[(wage['ROW'] == r) & (
                        wage['COL'] == c), race].item() * EPF * w_conv
                    
                    # Sort PM concentrations (Python is only lowest to highest :( )
                    pm_sort = np.sort(pm_out[str(Yr)][IC][(r, c)])
                    # So we reverse the order of the array elements
                    pm_descend = copy.deepcopy(pm_sort[::-1])
                    # Calculate dPM threshold for adaptation
                    dPM_thresh = (AT*w_i*10)/(base_mort*AF*VSL)
                    # Set adaptation count
                    count = 0
                    # Create deepcopy of pm concentrations in descending order
                    pm_test = copy.deepcopy(pm_descend)
    
                    # Calculate the potential reduction of mean outdoor PM by adapting each day
                    for i in range(len(pm_descend)):
                        pm_test[i] = ad_reduct * pm_test[i] 
                        pm_test_mean = np.mean(pm_test)
                        dPM_test = np.mean(pm_descend) - pm_test_mean
    
                        if dPM_test > dPM_thresh:
                            pm_descend[i] = pm_test[i]
                            count += 1
                        else:
                            pm_test[i] = pm_descend[i]
                            break
    
                    dPM = np.mean(pm_sort) - np.mean(pm_descend)
    
                    # Output adaptation days count
                    adapt_count[Yr][IC][(r, c)][race] = count
    
                    # Output the PM_out after adaptation
                    pm_out_adapt[Yr][IC][(r, c)][race] = np.mean(pm_descend)
    
                    # Output adaptation benefits
                    popu = pop.loc[(pop['ROW'] == r) & ( 
                        pop['COL'] == c), race].item() * PPF
                    adapt_ben[Yr][IC][(r, c)][race] = base_mort * AF * \
                        dPM * 0.1 * popu * ppa * VSL
                    adapt_ben_pc[Yr][IC][(r, c)][race] = base_mort * AF * dPM * 0.1 * VSL
    
                    # Output adaptation costs
                    adapt_hours = count * AT
                    adapt_cost[Yr][IC][(r, c)][race] = adapt_hours * w_i * popu * ppa
                    adapt_cost_pc[Yr][IC][(r, c)][race] = adapt_hours * w_i
    
                    # Output adaptation utility
                    adapt_util[Yr][IC][(r, c)][race] = adapt_ben[Yr][IC][(
                        r, c)][race] - adapt_cost[Yr][IC][(r, c)][race]
                    adapt_util_pc[Yr][IC][(r, c)][race] = adapt_ben_pc[Yr][IC][(
                        r, c)][race] - adapt_cost_pc[Yr][IC][(r, c)][race]


# Output pickles
# with open(outdir + 'adapt_ben_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_ben, f)
# f.close()

# with open(outdir + 'adapt_ben_pc_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_ben_pc, f)
# f.close()

# with open(outdir + 'adapt_cost_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_cost, f)
# f.close()

# with open(outdir + 'adapt_cost_pc_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_cost_pc, f)
# f.close()

# with open(outdir + 'adapt_count_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_count, f)
# f.close()

# with open(outdir + 'adapt_util_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_util, f)
# f.close()

# with open(outdir + 'adapt_util_pc_all_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(adapt_util_pc, f)
# f.close()

# with open(outdir + 'PM_out_base_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(pm_out_base, f)
# f.close()

# with open(outdir + 'PM_out_adapt_REF_' + str(inf_rate), 'wb') as f:
#     cp.dump(pm_out_adapt, f)
# f.close()
