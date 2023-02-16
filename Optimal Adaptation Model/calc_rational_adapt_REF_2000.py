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
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical Model PM\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"

# Load historical concentration data.
with open(indir2 + 'pm_conc_base_cell_model_2000', 'rb') as f:
     pm_out = cp.load(f)
f.close()

# Define list of years
years = list(range(1981, 2011))

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Import population per grid cell
pop = pd.read_csv(indir1 + "pop_lepeule.csv")

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Set adaptation time in hours. Assuming full adaptation because it maximizes utility.
AT = 2

#Set relative risk
RR = 1.14

# Calculate attributable fraction 
AF = (RR - 1)/ RR

# Set VSL
VSL = 7400000

# Calculate average annual exposure per cell
exp_all = {}
exp_mean = {}

for r, c in rc:
    exp_all[(r, c)] = []
    for Yr in years:
        exp_all[(r, c)].extend(pm_out[(r,c)][Yr])
        exp_mean[(r, c)] = np.mean(exp_all[(r, c)])

non_weighted_mean_conc = sum(exp_mean.values())/ 219

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
        pm_o = np.mean(pm_out[(r, c)][Yr])
        # Store baseline mean outdoor PM concentration in dictionary
        pm_out_base[(r, c)][Yr] = pm_o
        # Find and store mean hourly wage for the current grid cell
        w_i = wage.loc[(wage['ROW'] == r) & (wage['COL'] == c), '2000'].item()
        # Find and store baseline mortality rate for the current grid cell
        base_mort = bm.loc[(bm['ROW'] == r) & (bm['COL'] ==c), 'BaseMortRate'].item()
        # Set VSL
        VSL = 7400000
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
            pm_test[i] = 0
            pm_test_mean = np.mean(pm_test)
            dPM_test = np.mean(pm_descend) - pm_test_mean
            
            if dPM_test > dPM_thresh:
                pm_descend[i] = 0
                count += 1
            else:
                pm_test[i] = pm_descend[i]
                break
            
        dPM = np.mean(pm_sort) - np.mean(pm_descend)
        
        # Output adaptation days count    
        adapt_count[(r, c)][Yr] = count

        # Output the PM_out after adaptation
        pm_out_adapt[(r, c)][Yr] = np.mean(pm_descend)
        
        # Output adaptation benefits
        popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item()
        adapt_ben[(r, c)][Yr] = base_mort * AF * dPM * 0.1 * popu * VSL
        adapt_ben_pc[(r, c)][Yr] = base_mort * AF * dPM * 0.1 * VSL
        
        # Output adaptation costs
        adapt_hours = count * AT
        adapt_cost[(r, c)][Yr] = adapt_hours * w_i * popu 
        adapt_cost_pc[(r, c)][Yr] = adapt_hours * w_i  

        # Output adaptation utility
        adapt_util[(r, c)][Yr] = adapt_ben[(r, c)][Yr] - adapt_cost[(r, c)][Yr]
        adapt_util_pc[(r, c)][Yr] = adapt_ben_pc[(r, c)][Yr] - adapt_cost_pc[(r, c)][Yr]
                
# Calculate national level outcomes
# Pop-weighted adaptation exposure
pop_weighted_conc_ad = []
conc_dict = {}
for r, c in rc:
    conc_dict[(r, c)] = []
    for Yr in years:
        conc_dict[(r, c)].append(pm_out_adapt[(r, c)][Yr])
    conc_dict[(r, c)] = np.mean(conc_dict[(r, c)])
    pop_weighted_conc_ad.append(conc_dict[(r, c)] * pop_prop[(r, c)])

pop_weighted_conc_ad = np.sum(pop_weighted_conc_ad)

# Non pop weighted exposure
non_weighted_conc_ad = sum(conc_dict.values())/219

# Per capita base burden
pc_mort_base = pm_mort_nat_base / np.sum(pop['2000'])

# Total adaptation burden
pm_mort_val_ad = {}
for r, c in rc:
    bmr = bm.loc[(bm['ROW'] == r) & (bm['COL'] == c), 'BaseMortRate'].item()
    dPM = conc_dict[(r, c)]
    popu = pop.loc[(pop['ROW'] == r) * (pop['COL'] == c), '2000'].item()
    pm_mort_val_ad[(r, c)] = bmr * AF * dPM * 0.1 * popu * VSL

pm_mort_nat_ad = sum(pm_mort_val_ad.values())

# Per capita adaptation burden
pc_mort_nat_ad = pm_mort_nat_ad / np.sum(pop['2000'])

# Adaptation days per year
pop_weighted_ad_days = []
mean_ad_cell = {}
for r, c in rc:
    mean_ad_cell[(r, c)] = []
    for Yr in years:
        mean_ad_cell[(r, c)].append(adapt_count[(r, c)][Yr])
    mean_ad_cell[(r, c)] = np.mean(mean_ad_cell[(r, c)])
    pop_weighted_ad_days.append(mean_ad_cell[(r, c)] * pop_prop[(r, c)])

pop_weighted_ad_days = round(np.sum(pop_weighted_ad_days))

# Calculate national mean total adaptation utility
nat_ad_util = []
for r, c in rc:
    nat_ad_util.append(sum(adapt_util[(r, c)].values())/30)

nat_ad_util = np.sum(nat_ad_util)/ 1E12

                
# Output pickles 
# with open(outdir + 'adapt_ben_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_ben,f)
# f.close()

# with open(outdir + 'adapt_ben_pc_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_ben_pc,f)
# f.close()

# with open(outdir + 'adapt_cost_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_cost,f)
# f.close()

# with open(outdir + 'adapt_cost_pc_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_cost_pc,f)
# f.close()

# with open(outdir + 'adapt_count_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_count,f)
# f.close()

# with open(outdir + 'adapt_util_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_util,f)
# f.close()

# with open(outdir + 'adapt_util_pc_all_REF_2000', 'wb') as f:
#       cp.dump(adapt_util_pc,f)
# f.close()

# with open(outdir + 'PM_out_adapt_REF_2000', 'wb') as f:
#       cp.dump(pm_out_adapt,f)
# f.close()

# with open(outdir + 'PM_out_base_REF_2000', 'wb') as f:
#       cp.dump(pm_out_base,f)
# f.close()
