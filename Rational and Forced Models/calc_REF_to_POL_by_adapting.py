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

inf_rate = 0.2

# Set directories
indir1 = r".\Data\General\\"
indir2 = r".\Data\Future AQ\\"
indir3 = r".\Data\AGU2020-Fig3\\"
outdir = r".\Data\Rational Actor Infiltration\\"

# Load future concentration data. REF, P3.7 or P4.5
with open(indir2 + 'REF_PM25_2050_2100_d', 'rb') as f:
    pm_out_REF = cp.load(f)
f.close()

with open(indir2 + 'P45_PM25_2050_2100_d', 'rb') as f:
    pm_out_P45 = cp.load(f)
f.close()

with open(indir2 + 'P37_PM25_2050_2100_d', 'rb') as f:
    pm_out_P37 = cp.load(f)
f.close()


# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis=0, inplace=True)
# Reset index back to starting at 0
wage.reset_index(inplace=True)

# Import population per grid cell
pop = pd.read_csv(indir1 + "pop_lepeule.csv")

# Create list of all grid cells
rc = list(zip(list(wage['ROW']), list(wage['COL'])))

# Calculate proportion of population in each grid cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')


# Set adaptation time in hours. Assuming full adaptation because it maximizes utility.
AT = 2

# Set relative risk
RR = 1.14

# Calculate attributable fraction
AF = (RR - 1) / RR

# Set VSL
VSL = 7400000

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

# Calculate adaptation days for all scenarios
# Initialize dictionaries
pm_out_base = {}
pm_out_pol = {}
adapt_count = {}


# Iterate over both sets of years (for 2050 and 2100)
for pol in ['P37','P45']:
    adapt_count[pol] = {}
    for years in year_set:
        # Iterate over all years in year set
        if np.mean(years) < 2070:
            # Economic projection factor
            EPF = 2.01
            # Population projection factor
            PPF = 1.48
            # PM Mortality Factor (for baseline mortaltiy)
            PMMF = 1.12
            # Set VSL
            VSL = 7400000
            # Update VSL based on EPF
            VSL = VSL * EPF
        else:
            EPF = 3.94
            PPF = 1.73
            PMMF = 0.92
            # Set VSL
            VSL = 7400000
            # Update VSL based on EPF
            VSL = VSL * EPF
        for Yr in years:
            # Initialize next level of dictionaries
            pm_out_base[Yr] = {}
            pm_out_pol[Yr] = {}
            adapt_count[pol][Yr] = {}
            
    
            # Iterate over all initial conditions in each year
            for IC in IC_list:
                # Initialize next level of dictionaries
                pm_out_base[Yr][IC] = {}
                pm_out_pol[Yr][IC] = {}
                adapt_count[pol][Yr][IC] = {}
    
                # Iterate over all cells in each scenario
                for r, c in rc:
                    # Calculate pol mean outdoor PM concentration
                    if pol == 'P37':
                        pm_target = np.mean(pm_out_P37[str(Yr)][IC][(r, c)])
                    elif pol == 'P45':
                        pm_target = np.mean(pm_out_P45[str(Yr)][IC][(r, c)])
                    
                    
                    # Calculate baseline mean outdoor PM concentration
                    pm_o = np.mean(pm_out_REF[str(Yr)][IC][(r, c)])
                    # Store baseline mean outdoor PM concentration in dictionary
                    pm_out_base[Yr][IC][(r, c)] = pm_o
                    # Sort PM concentrations (Python is only lowest to highest :( )
                    pm_sort = np.sort(pm_out_REF[str(Yr)][IC][(r, c)])
                    # So we reverse the order of the array elements
                    pm_descend = copy.deepcopy(pm_sort[::-1])
                    # Calculate dPM threshold for adaptation
                    dPM_thresh = pm_o - pm_target
                    # Set adaptation count
                    count = 0
                    # Create deepcopy of pm concentrations in descending order
                    pm_test = copy.deepcopy(pm_descend)
                    # Creat dPM_total
                    dPM_total = 0
    
                    # Calculate the potential reduction of mean outdoor PM by adapting each day
                    for i in range(len(pm_descend)):
                        pm_test[i] = ad_reduct * pm_test[i]
                        pm_test_mean = np.mean(pm_test)
                        dPM_test = np.mean(pm_descend) - pm_test_mean
                        
                        
                        if dPM_thresh <= 0:
                            break
                        elif dPM_total < dPM_thresh:
                            pm_descend[i] = pm_test[i]
                            dPM_total += dPM_test
                            count += 1
                        else:
                            break
    
                    # Output adaptation days count
                    adapt_count[pol][Yr][IC][(r, c)] = count
    
#%%
# Convert adapt count by cell to national pop-weighted average
adapt_nat = {}


adapt_cell_2050 = {}
adapt_cell_2100 = {}

for pol in ['P37', 'P45']:
    adapt_cell_2050[pol] = {}
    adapt_cell_2100[pol] = {}
    for r, c in rc:
        adapt_cell_2050[pol][(r, c)] = []
        adapt_cell_2100[pol][(r, c)] = []
        for years in year_set:
            for Yr in years:
                for IC in IC_list:
                    if Yr < 2070:
                        adapt_cell_2050[pol][(r, c)].append(adapt_count[pol][Yr][IC][(r, c)])
                    else:
                        adapt_cell_2100[pol][(r, c)].append(adapt_count[pol][Yr][IC][(r, c)])
                        
# Calculate range of national days
for year in [2050, 2100]:
    adapt_nat[year] = {}
    for pol in ['P37', 'P45']:
        adapt_nat[year][pol] = []
        if year == 2050:
            for i in range(0, 150):
                year_ad = []
                for r, c in rc:
                    prop_ref = pop_prop[(r, c)]
                    ad_cell = adapt_cell_2050[pol][(r, c)][i]
                    year_ad.append(prop_ref * ad_cell)
                adapt_nat[year][pol].append(round(sum(year_ad)))    
        else:
            for i in range(0, 150):
                year_ad = []
                for r, c in rc:
                    prop_ref = pop_prop[(r, c)]
                    ad_cell = adapt_cell_2100[pol][(r, c)][i]
                    year_ad.append(prop_ref * ad_cell)
                adapt_nat[year][pol].append(round(sum(year_ad))) 
            
mean_adapt_nat = {}
for year in [2050,2100]:
    mean_adapt_nat[year] = {}
    for pol in ['P37','P45']:
        mean_adapt_nat[year][pol] = round(np.mean(adapt_nat[year][pol]))
       


# Convert to 2020 dollars # Note not factored to future economies
mean_wage = 14.11 
# Calculate national
cost_adapt_nat = {}
for year in [2050,2100]:
    cost_adapt_nat[year] = {}
    for pol in ['P37','P45']:
        cost_adapt_nat[year][pol] = mean_wage * mean_adapt_nat[year][pol]

# Output pickle