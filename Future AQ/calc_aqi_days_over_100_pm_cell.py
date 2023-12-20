# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 09:01:45 2022

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
indir1 = r".\Data\Future AQ\\"
indir2 = r".\Data\BenMap Grid\\"
outdir = r".\Data\Future AQ\\"

# Load Future AQI data
with open(indir1 + 'REF_prop_bad_2050_2100_d', 'rb') as f:
     prop_REF = cp.load(f)
f.close()

with open(indir1 + 'P45_prop_bad_2050_2100_d', 'rb') as f:
     prop_P45 = cp.load(f)
f.close()

with open(indir1 + 'P37_prop_bad_2050_2100_d', 'rb') as f:
     prop_P37 = cp.load(f)
f.close()

# Load BenMap Grid
grid = gpd.read_file(indir2 + "BenMapGridPoints.csv")
grid = grid.drop('geometry', axis = 1)
grid['geom'] = grid['geom'].apply(wkt.loads)
grid = grid.set_geometry('geom')
grid['ROW'] = grid['ROW'].astype('int')
grid['COL'] = grid['COL'].astype('int')
grid = grid.drop([1, 219, 220], axis = 0)

grid2050 = copy.deepcopy(grid)
grid2100 = copy.deepcopy(grid)
# Create rc
rc = list(zip(grid['ROW'], grid['COL']))

# Fill in grid with 150 year means 2050
grid2050['REF'] = 0
grid2050['P45'] = 0
grid2050['P37'] = 0

for r, c in rc:
    prop_REF_list = []
    prop_P37_list = []
    prop_P45_list = []
    for Yr in range(2036, 2066):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            prop_REF_list.append(prop_REF[str(Yr)][IC][(r, c)])
            prop_P37_list.append(prop_P37[str(Yr)][IC][(r, c)])
            prop_P45_list.append(prop_P45[str(Yr)][IC][(r, c)])
            
    grid2050.loc[(grid2050['ROW'] == r) & (grid2050['COL'] == c), 'REF'] = np.mean(prop_REF_list)
    grid2050.loc[(grid2050['ROW'] == r) & (grid2050['COL'] == c), 'P45'] = np.mean(prop_P45_list)
    grid2050.loc[(grid2050['ROW'] == r) & (grid2050['COL'] == c), 'P37'] = np.mean(prop_P37_list)

# Fill in grid with 150 year means 2100
grid2100['REF'] = 0
grid2100['P45'] = 0
grid2100['P37'] = 0

for r, c in rc:
    prop_REF_list = []
    prop_P37_list = []
    prop_P45_list = []
    for Yr in range(2086, 2116):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            prop_REF_list.append(prop_REF[str(Yr)][IC][(r, c)])
            prop_P37_list.append(prop_P37[str(Yr)][IC][(r, c)])
            prop_P45_list.append(prop_P45[str(Yr)][IC][(r, c)])
            
    grid2100.loc[(grid2100['ROW'] == r) & (grid2100['COL'] == c), 'REF'] = np.mean(prop_REF_list)
    grid2100.loc[(grid2100['ROW'] == r) & (grid2100['COL'] == c), 'P45'] = np.mean(prop_P45_list)
    grid2100.loc[(grid2100['ROW'] == r) & (grid2100['COL'] == c), 'P37'] = np.mean(prop_P37_list)


grid2050['REF Days'] = round(grid2050['REF'] * 365)
grid2050['P45 Days'] = round(grid2050['P45'] * 365)
grid2050['P37 Days'] = round(grid2050['P37'] * 365)

grid2100['REF Days'] = round(grid2100['REF'] * 365)
grid2100['P45 Days'] = round(grid2100['P45'] * 365)
grid2100['P37 Days'] = round(grid2100['P37'] * 365)


# Output pickles of grid2050 and 2100
with open(outdir + 'grid_days_bad_2050_pm_100', 'wb') as f:
      cp.dump(grid2050,f)
f.close()

with open(outdir + 'grid_days_bad_2100_pm_100', 'wb') as f:
      cp.dump(grid2100,f)
f.close()