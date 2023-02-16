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
with open(indir1 + 'adapt_ben_all_P37', 'rb') as f:
     ben_P37 = cp.load(f)
f.close()

with open(indir1 + 'adapt_ben_all_P45', 'rb') as f:
     ben_P45 = cp.load(f)
f.close()

with open(indir1 + 'adapt_ben_all_REF', 'rb') as f:
     ben_REF = cp.load(f)
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
    
# Calculate national costs 2050 by cell
ben_cell_REF_2050 = {}
ben_cell_P45_2050 = {}
ben_cell_P37_2050 = {}

for r, c in rc:
    ben_cell_REF_2050[(r, c)] = []
    ben_cell_P45_2050[(r, c)] = []
    ben_cell_P37_2050[(r, c)] = []
    
for Yr in year_set[0]:
    for IC in IC_list:
        for r, c in rc:
            ben_cell_REF_2050[(r, c)].append(ben_REF[Yr][IC][(r, c)])
            ben_cell_P45_2050[(r, c)].append(ben_P45[Yr][IC][(r, c)])
            ben_cell_P37_2050[(r, c)].append(ben_P37[Yr][IC][(r, c)])

for r, c in rc:            
    ben_cell_REF_2050[(r, c)] = np.mean(ben_cell_REF_2050[(r, c)])
    ben_cell_P45_2050[(r, c)] = np.mean(ben_cell_P45_2050[(r, c)])
    ben_cell_P37_2050[(r, c)] = np.mean(ben_cell_P37_2050[(r, c)])
    
nat_ben_REF_2050 = sum(ben_cell_REF_2050.values()) / 1E12
nat_ben_P45_2050 = sum(ben_cell_P45_2050.values()) / 1E12
nat_ben_P37_2050 = sum(ben_cell_P37_2050.values()) / 1E12

# Calculate national costs 2100 by cell
ben_cell_REF_2100 = {}
ben_cell_P45_2100 = {}
ben_cell_P37_2100 = {}

for r, c in rc:
    ben_cell_REF_2100[(r, c)] = []
    ben_cell_P45_2100[(r, c)] = []
    ben_cell_P37_2100[(r, c)] = []
    
for Yr in year_set[1]:
    for IC in IC_list:
        for r, c in rc:
            ben_cell_REF_2100[(r, c)].append(ben_REF[Yr][IC][(r, c)])
            ben_cell_P45_2100[(r, c)].append(ben_P45[Yr][IC][(r, c)])
            ben_cell_P37_2100[(r, c)].append(ben_P37[Yr][IC][(r, c)])

for r, c in rc:            
    ben_cell_REF_2100[(r, c)] = np.mean(ben_cell_REF_2100[(r, c)])
    ben_cell_P45_2100[(r, c)] = np.mean(ben_cell_P45_2100[(r, c)])
    ben_cell_P37_2100[(r, c)] = np.mean(ben_cell_P37_2100[(r, c)])
    
nat_ben_REF_2100 = sum(ben_cell_REF_2100.values()) / 1E12
nat_ben_P45_2100 = sum(ben_cell_P45_2100.values()) / 1E12
nat_ben_P37_2100 = sum(ben_cell_P37_2100.values()) / 1E12


