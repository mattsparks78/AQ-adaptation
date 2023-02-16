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
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"

# Load pickles
# 2050
with open(indir1 + 'pm_base_REF_cell_2050', 'rb') as f:
     pm_b_REF_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2050', 'rb') as f:
     pm_ad_REF_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P45_cell_2050', 'rb') as f:
     pm_b_P45_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P45_cell_2050', 'rb') as f:
     pm_ad_P45_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2050', 'rb') as f:
     pm_b_P37_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2050', 'rb') as f:
     pm_ad_P37_2050 = cp.load(f)
f.close()

# 2100
with open(indir1 + 'pm_base_REF_cell_2100', 'rb') as f:
     pm_b_REF_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2100', 'rb') as f:
     pm_ad_REF_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P45_cell_2100', 'rb') as f:
     pm_b_P45_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P45_cell_2100', 'rb') as f:
     pm_ad_P45_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2100', 'rb') as f:
     pm_b_P37_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2100', 'rb') as f:
     pm_ad_P37_2100 = cp.load(f)
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
nat_REF_2050_ad = []
nat_P45_2050_b = []
nat_P45_2050_ad = []
nat_P37_2050_b = []
nat_P37_2050_ad = []

nat_REF_2100_b = []
nat_REF_2100_ad = []
nat_P45_2100_b = []
nat_P45_2100_ad = []
nat_P37_2100_b = []
nat_P37_2100_ad = []

for r, c in rc:
    nat_REF_2050_b.append(np.mean(pm_b_REF_2050[(r, c)]) * pop_prop[(r, c)])
    nat_REF_2050_ad.append(np.mean(pm_ad_REF_2050[(r, c)]) * pop_prop[(r, c)])
    nat_P45_2050_b.append(np.mean(pm_b_P45_2050[(r, c)]) * pop_prop[(r, c)])
    nat_P45_2050_ad.append(np.mean(pm_ad_P45_2050[(r, c)]) * pop_prop[(r, c)])
    nat_P37_2050_b.append(np.mean(pm_b_P37_2050[(r, c)]) * pop_prop[(r, c)])
    nat_P37_2050_ad.append(np.mean(pm_ad_P37_2050[(r, c)]) * pop_prop[(r, c)])

    nat_REF_2100_b.append(np.mean(pm_b_REF_2100[(r, c)]) * pop_prop[(r, c)])
    nat_REF_2100_ad.append(np.mean(pm_ad_REF_2100[(r, c)]) * pop_prop[(r, c)])
    nat_P45_2100_b.append(np.mean(pm_b_P45_2100[(r, c)]) * pop_prop[(r, c)])
    nat_P45_2100_ad.append(np.mean(pm_ad_P45_2100[(r, c)]) * pop_prop[(r, c)])
    nat_P37_2100_b.append(np.mean(pm_b_P37_2100[(r, c)]) * pop_prop[(r, c)])
    nat_P37_2100_ad.append(np.mean(pm_ad_P37_2100[(r, c)]) * pop_prop[(r, c)])

nat_REF_2050_b = sum(nat_REF_2050_b)
nat_REF_2050_ad = sum(nat_REF_2050_ad)
nat_P45_2050_b = sum(nat_P45_2050_b)
nat_P45_2050_ad = sum(nat_P45_2050_ad)
nat_P37_2050_b = sum(nat_P37_2050_b)
nat_P37_2050_ad = sum(nat_P37_2050_ad)

nat_REF_2100_b = sum(nat_REF_2100_b)
nat_REF_2100_ad = sum(nat_REF_2100_ad)
nat_P45_2100_b = sum(nat_P45_2100_b)
nat_P45_2100_ad = sum(nat_P45_2100_ad)
nat_P37_2100_b = sum(nat_P37_2100_b)
nat_P37_2100_ad = sum(nat_P37_2100_ad)
