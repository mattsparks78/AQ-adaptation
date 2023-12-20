# -*- coding: utf-8 -*-
"""
Created on Sun May 29 15:41:29 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp



# Load total SCI dictionary
with open(r".\Data\SCI\sci_tot_b", 'rb') as f:
     sci_tot = cp.load(f)
f.close()

# Load population per celldata
grid_pop = pd.read_csv(r".\Data\BenMap Grid\BenMapGridPop.csv")
grid_pop = grid_pop[['ROW', 'COL', '2000']]
grid_pop.drop([1, 219, 220], axis = 0, inplace = True)

# Create list of cells
rc = list(zip(list(grid_pop['ROW']), list(grid_pop['COL'])))

# Create dictionary with population values
pop_dict = {}
for r, c in rc:
    pop_dict[(r, c)] = grid_pop.loc[(grid_pop['ROW'] == r) & (grid_pop['COL'] == c), '2000'].item()

# Calculate connections = SCI * pop_i * pop_j
con = {}
for r, c in rc:
    con[(r, c)] = {}
    for row, col in rc:
        con[(r, c)][(row, col)] = sci_tot[(r, c)][(row, col)] * pop_dict[(r, c)] * pop_dict[(row, col)]

# Normalize cell values by total number of connections
sci_norm = {}
for r, c in rc:
    sci_norm[(r,c)] = {}
    for row, col in rc:
        sci_norm[(r, c)][(row, col)] = con[(r, c)][(row, col)] / sum(con[(r, c)].values())


# Save normalized SCI dictionary      
with open(r'.\Data\SCI\sci_normalized_b', 'wb') as f:
     cp.dump(sci_norm,f)
f.close()
        