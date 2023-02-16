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
with open(indir1 + 'mean_adapt_cell_REF_2050', 'rb') as f:
     ad_REF_2050 = cp.load(f)
f.close()

with open(indir1 + 'mean_adapt_cell_P45_2050', 'rb') as f:
     ad_P45_2050 = cp.load(f)
f.close()

with open(indir1 + 'mean_adapt_cell_P37_2050', 'rb') as f:
     ad_P37_2050 = cp.load(f)
f.close()

# 2100
with open(indir1 + 'mean_adapt_cell_REF_2100', 'rb') as f:
     ad_REF_2100 = cp.load(f)
f.close()

with open(indir1 + 'mean_adapt_cell_P45_2100', 'rb') as f:
     ad_P45_2100 = cp.load(f)
f.close()

with open(indir1 + 'mean_adapt_cell_P37_2100', 'rb') as f:
     ad_P37_2100 = cp.load(f)
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
    
# Calculate national population weighted adaptation days for each scenario
ad_days_REF_2050 = []
ad_days_P45_2050 = []
ad_days_P37_2050 = []
ad_days_REF_2100 = []
ad_days_P45_2100 = []
ad_days_P37_2100 = []

for r, c in rc:
    ad_days_REF_2050.append(ad_REF_2050[(r, c)] * pop_prop[(r, c)])
    ad_days_P45_2050.append(ad_P45_2050[(r, c)] * pop_prop[(r, c)])
    ad_days_P37_2050.append(ad_P37_2050[(r, c)] * pop_prop[(r, c)])
    ad_days_REF_2100.append(ad_REF_2100[(r, c)] * pop_prop[(r, c)])
    ad_days_P45_2100.append(ad_P45_2100[(r, c)] * pop_prop[(r, c)])
    ad_days_P37_2100.append(ad_P37_2100[(r, c)] * pop_prop[(r, c)])

ad_days_REF_2050 = round(np.sum(ad_days_REF_2050))
ad_days_P45_2050 = round(np.sum(ad_days_P45_2050))
ad_days_P37_2050 = round(np.sum(ad_days_P37_2050))
ad_days_REF_2100 = round(np.sum(ad_days_REF_2100))
ad_days_P45_2100 = round(np.sum(ad_days_P45_2100))
ad_days_P37_2100 = round(np.sum(ad_days_P37_2100))

days_dict = {}
days_dict['REF'] = {}
days_dict['P45'] = {}
days_dict['P37'] = {}

days_dict['REF'][2050] = ad_days_REF_2050
days_dict['REF'][2100] = ad_days_REF_2100
days_dict['P45'][2050] = ad_days_P45_2050
days_dict['P45'][2100] = ad_days_P45_2100
days_dict['P37'][2050] = ad_days_P37_2050
days_dict['P37'][2100] = ad_days_P37_2100





