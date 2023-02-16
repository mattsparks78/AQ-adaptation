# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 09:47:14 2022

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

# Load historical data.
with open(indir1 + 'adapt_cost_all_REF_2000', 'rb') as f:
     cost_all = cp.load(f)
f.close()

with open(indir1 + 'adapt_ben_all_REF_2000', 'rb') as f:
     ben_all = cp.load(f)
f.close()

with open(indir1 + 'adapt_util_all_REF_2000', 'rb') as f:
     util_all = cp.load(f)
f.close()

# Define list of years
years = list(range(1981, 2011))

# Import hourly wage data
wage = pd.read_csv(indir2 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Import population per grid cell
pop = pd.read_csv(indir2 + "pop_lepeule.csv")

# Import baseline mortality per grid cell
bm = pd.read_csv(indir2 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Aggregate cells to national level
nat_ben = []
nat_cost = []
nat_util = []

for r, c in rc:
    nat_ben.append(sum(ben_all[(r, c)].values()) / 30)
    nat_cost.append(sum(cost_all[(r, c)].values()) / 30)
    nat_util.append(sum(util_all[(r, c)].values()) / 30)

nat_ben = np.sum(nat_ben) / 1E12
nat_cost = np.sum(nat_cost) / 1E12
nat_util = np.sum(nat_util) / 1E12
                    

