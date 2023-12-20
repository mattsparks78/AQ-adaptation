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
indir2 = r".\Data\Historical Model PM\\"
indir3 = r".\Data\AGU2020-Fig3\\"
outdir = r".\Data\Rational Actor Infiltration Dem\\"

# Load historical concentration data.
with open(indir2 + 'pm_conc_base_cell_model_2000', 'rb') as f:
     pm_out = cp.load(f)
f.close()

# Define list of years
years = list(range(1981, 2011))

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_race.csv")

# Import population data
pop = pd.read_csv(indir1 + 'pop_lepeule_race.csv')

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Calculate fractions of indoor and outdoor time
t_out = 1
t_in = 24 - t_out

prop_out = t_out / 24
prop_in = t_in / 24

# Set adaptation time in hours. 
AT = 1

# Set infiltration rate based on empirical data/for sensitvity analysis
inf_rate = 0.2

# Calculate ratio of adaptation exposure to baseline exposure
exp_base = prop_out + prop_in*inf_rate
exp_ad = (prop_out + prop_in)*inf_rate

ad_reduct = exp_ad / exp_base
#Set relative risk
RR = 1.14

# Calculate attributable fraction 
AF = (RR - 1)/ RR

# Set VSL. Converted from $2006 USD to 2005
VSL = 7104000

# Set conversion factor for 2000 to 2005 USD for wage
w_conv = 1.13

# Calculate average annual exposure per cell
exp_all = {}
exp_mean = {}

for r, c in rc:
    exp_all[(r, c)] = []
    for Yr in years:
        exp_all[(r, c)].extend(pm_out[(r,c)][Yr])
        exp_mean[(r, c)] = np.mean(exp_all[(r, c)])

# non_weighted_mean_conc = sum(exp_mean.values())/ 219

# Calculate proportion of population in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    
# Calculate population weighted national exposure
pop_weighted_conc = []
for r, c in rc:
    pop_weighted_conc.append(exp_mean[(r, c)] * pop_prop[(r, c)])

pop_weighted_conc = np.sum(pop_weighted_conc)

# Calculate base health burden from PM
pm_mort_val_base = {}
for r, c in rc:
    bmr = bm.loc[(bm['ROW'] == r) & (bm['COL'] == c), 'BaseMortRate'].item()
    dPM = exp_mean[(r, c)]
    popu = pop.loc[(pop['ROW'] == r) * (pop['COL'] == c), '2000'].item()
    pm_mort_val_base[(r, c)] = bmr * AF * dPM * 0.1 * popu * VSL

# Calculate national PM mort pop weighted
pm_mort_nat_base = sum(pm_mort_val_base.values())

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
for r, c in rc:
    # Initialize next level of dictionaries
    pm_out_base[(r, c)] = {}
    pm_out_adapt[(r, c)] = {}
    adapt_count[(r, c)] = {}
    adapt_cost[(r, c)] = {}
    adapt_ben[(r, c)] = {}
    adapt_util[(r, c)] = {}
    adapt_cost_pc[(r, c)] = {}
    adapt_ben_pc[(r, c)] = {}
    adapt_util_pc[(r, c)] = {}
        
    # Iterate over all years
    for Yr in years:
        pm_out_adapt[(r, c)][Yr] = {}
        adapt_count[(r, c)][Yr] = {}
        adapt_cost[(r, c)][Yr] = {}
        adapt_ben[(r, c)][Yr] = {}
        adapt_util[(r, c)][Yr] = {}
        adapt_cost_pc[(r, c)][Yr] = {}
        adapt_ben_pc[(r, c)][Yr] = {}
        adapt_util_pc[(r, c)][Yr] = {}
        
        pm_o = np.mean(pm_out[(r, c)][Yr])
        # Store baseline mean outdoor PM concentration in dictionary
        pm_out_base[(r, c)][Yr] = pm_o
        # Find and store baseline mortality rate for the current grid cell
        base_mort = bm.loc[(bm['ROW'] == r) & (bm['COL'] ==c), 'BaseMortRate'].item()
        
        # Iterate over all races
        for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
            # Find and store mean hourly wage for the current grid cell
            w_i = wage.loc[(wage['ROW'] == r) & (wage['COL'] == c), race].item() * w_conv
            
        
       
            # Sort PM concentrations (Python is only lowest to highest :( )
            pm_sort = np.sort(pm_out[(r, c)][Yr])
            # So we reverse the order of the array elements
            pm_descend = copy.deepcopy(pm_sort[::-1])
            # Calculate total PM
            pm_tot = np.sum(pm_descend)
            # Calculate dPM threshold for adaptation
            dPM_thresh = (AT*w_i*10)/(base_mort*AF*VSL)
            # Set adaptation count
            count = 0
        
        
            # Calculate the potential reduction of mean outdoor PM by adapting each day
            for i in range(len(pm_descend)):
                pm_test = copy.deepcopy(pm_descend)
                pm_test[i] = ad_reduct * pm_test[1]
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
            adapt_count[(r, c)][Yr][race] = count

            # Output the PM_out after adaptation
            pm_out_adapt[(r, c)][Yr][race] = np.mean(pm_descend)
        
            # Output adaptation benefits
            popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), race].item()
            adapt_ben[(r, c)][Yr][race] = base_mort * AF * dPM * 0.1 * popu * VSL
            adapt_ben_pc[(r, c)][Yr][race] = base_mort * AF * dPM * 0.1 * VSL
        
            # Output adaptation costs
            adapt_hours = count * AT
            adapt_cost[(r, c)][Yr][race] = adapt_hours * w_i * popu 
            adapt_cost_pc[(r, c)][Yr][race] = adapt_hours * w_i  

            # Output adaptation utility
            adapt_util[(r, c)][Yr][race] = adapt_ben[(r, c)][Yr][race] - adapt_cost[(r, c)][Yr][race]
            adapt_util_pc[(r, c)][Yr][race] = adapt_ben_pc[(r, c)][Yr][race] - adapt_cost_pc[(r, c)][Yr][race]

                
# # Output pickles 
# with open(outdir + 'adapt_ben_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_ben,f)
# f.close()

# with open(outdir + 'adapt_ben_pc_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_ben_pc,f)
# f.close()

# with open(outdir + 'adapt_cost_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_cost,f)
# f.close()

# with open(outdir + 'adapt_cost_pc_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_cost_pc,f)
# f.close()

# with open(outdir + 'adapt_count_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_count,f)
# f.close()

# with open(outdir + 'adapt_util_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_util,f)
# f.close()

# with open(outdir + 'adapt_util_pc_all_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_util_pc,f)
# f.close()

# with open(outdir + 'PM_out_adapt_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_out_adapt,f)
# f.close()

# with open(outdir + 'PM_out_base_REF_2000_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_out_base,f)
# f.close()
