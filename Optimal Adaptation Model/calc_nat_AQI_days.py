# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 09:21:37 2022

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
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical Model PM\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"

# Load pickles
# 2050
with open(indir1 + 'REF_tot_bad_2050_2100_d_pm_only', 'rb') as f:
     pm_b_REF = cp.load(f)
f.close()

with open(indir1 + 'P45_tot_bad_2050_2100_d_pm_only', 'rb') as f:
     pm_b_P45 = cp.load(f)
f.close()

with open(indir1 + 'P37_tot_bad_2050_2100_d_pm_only', 'rb') as f:
     pm_b_P37 = cp.load(f)
f.close()

# 2000
with open(indir3 + 'pm_conc_base_cell_model_2000', 'rb') as f:
     pm_2000 = cp.load(f)
f.close()

# Load remaining things
# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir2 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Import population per grid cell
pop = pd.read_csv(indir2 + "pop_lepeule.csv")

# Calculate proportion of population in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    
# Calculate national population weighted average concentration for each scenario
nat_REF_2050_b = []
nat_P45_2050_b = []
nat_P37_2050_b = []

nat_REF_2100_b = []
nat_P45_2100_b = []
nat_P37_2100_b = []

for year in year_set[0]:
    for IC in IC_list:
        IC_list_REF = []
        IC_list_P45 = []
        IC_list_P37 = []
        for r, c in rc:
            IC_list_REF.append(pm_b_REF[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            IC_list_P45.append(pm_b_P45[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            IC_list_P37.append(pm_b_P37[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            
        nat_REF_2050_b.append(round(sum(IC_list_REF)))
        nat_P45_2050_b.append(round(sum(IC_list_P45)))
        nat_P37_2050_b.append(round(sum(IC_list_P37)))

REF_2050_min = min(nat_REF_2050_b)
REF_2050_max = max(nat_REF_2050_b)

P45_2050_min = min(nat_P45_2050_b)
P45_2050_max = max(nat_P45_2050_b)

P37_2050_min = min(nat_P37_2050_b)
P37_2050_max = max(nat_P37_2050_b)
    
for year in year_set[1]:
    for IC in IC_list:
        IC_list_REF = []
        IC_list_P45 = []
        IC_list_P37 = []
        for r, c in rc:
            IC_list_REF.append(pm_b_REF[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            IC_list_P45.append(pm_b_P45[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            IC_list_P37.append(pm_b_P37[str(year)][IC][(r, c)] * pop_prop[(r, c)])
            
        nat_REF_2100_b.append(round(sum(IC_list_REF)))
        nat_P45_2100_b.append(round(sum(IC_list_P45)))
        nat_P37_2100_b.append(round(sum(IC_list_P37)))

REF_2100_min = min(nat_REF_2100_b)
REF_2100_max = max(nat_REF_2100_b)

P45_2100_min = min(nat_P45_2100_b)
P45_2100_max = max(nat_P45_2100_b)

P37_2100_min = min(nat_P37_2100_b)
P37_2100_max = max(nat_P37_2100_b)


# Calc 2000. Calc bad air days 
bad_2000 = {}

for r, c in rc:
    bad_2000[(r, c)] = {}
    for year in range(1981, 2011):
        count_bad = len([x for x in pm_2000[(r, c)][year] if x > 35.5])
        bad_2000[(r, c)][year] = count_bad
# Calc national annual values max/min
nat_list_2000 = []
for year in range(1981, 2011):
    year_list = []
    for r,c in rc:
        year_list.append(bad_2000[(r, c)][year] * pop_prop[(r, c)])
    nat_list_2000.append(round(sum(year_list)))