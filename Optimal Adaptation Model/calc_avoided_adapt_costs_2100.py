# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 11:30:08 2022

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

# Load necessary data
with open(indir1 + 'mean_adapt_cell_P37_2100', 'rb') as f:
    adapt_P37 = cp.load(f)
f.close()

with open(indir1 + 'mean_adapt_cell_REF_2100', 'rb') as f:
    adapt_REF = cp.load(f)
f.close()

# Calculate difference in adaptation days
# Load wage data
wage = pd.read_csv(indir2 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis=0, inplace=True)
# Reset index back to starting at 0
wage.reset_index(inplace=True)

# Create list of all grid cells
rc = list(zip(list(wage['ROW']), list(wage['COL'])))

# Import population per grid cell
pop = pd.read_csv(indir2 + "benmap_pop_all.csv")
# Drop grid cells that aren't included in BenMap analysis
pop.drop([1, 219, 220], axis=0, inplace=True)
# Reset index back to starting at 0
pop.reset_index(inplace=True)

# Economic and Population Factors
EPF = 3.37
PPF = 1.73

dif_ad = {}
for r, c in rc:
    dif_ad[(r, c)] = adapt_REF[(r, c)] - adapt_P37[(r, c)]

dif_cost = {}
for r, c in rc:
    popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * EPF
    wages = wage.loc[(wage['ROW'] == r) & (
        wage['COL'] == c), '2000'].item() * PPF
    dif_cost[(r, c)] = dif_ad[(r, c)] * 1.92 * popu * wages

ad_cost_dif = sum(dif_cost.values())
ad_cost_dif / 1E9




