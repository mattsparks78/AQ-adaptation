# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:31:36 2022

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
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical Model\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical Model\\"

# Define parameters
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
pop = pd.read_csv(indir1 + "benmap_pop_all.csv")
# Drop grid cells that aren't included in BenMap analysis
pop.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
pop.reset_index(inplace = True)

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Calculate proportion of population in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    
# Load PM25 data and place in dictionary
# Initialize dictionary
pm_conc_base = {}

# Add data to dictionary
for Yr in years:
    pm_conc_base[Yr] = {}
    filename = 'REF.CS30.MIC.1981-2010.' + str(Yr) + '.AQ.csv'
    pm_data = pd.read_csv(indir2 + filename, usecols = ['Column', 'Row', 'Values'])
    for r,c in rc:
        pm_conc_base[Yr][(r, c)] = pm_data.loc[(pm_data['Row'] == r) & (pm_data['Column'] == c), 'Values'].item()
        pm_conc_base[Yr][(r, c)] = pm_conc_base[Yr][(r, c)].split(',')
        pm_conc_base[Yr][(r, c)] = [float(x) for x in pm_conc_base[Yr][(r, c)]]

# Break down data by cell
pm_conc_base_cell = {}

# Add data to dictionary
for r, c in rc:
    pm_conc_base_cell[(r, c)] = {}
    for Yr in years:
        pm_conc_base_cell[(r, c)][Yr] = pm_conc_base[Yr][(r, c)]
        
# Output pickles
with open(outdir + "pm_conc_base_model_2000", 'wb') as f:
      cp.dump(pm_conc_base,f)
f.close()

with open(outdir + "pm_conc_base_cell_model_2000", 'wb') as f:
      cp.dump(pm_conc_base_cell,f)
f.close()
