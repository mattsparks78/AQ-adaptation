# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 14:31:17 2022

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
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"

# Load per capita adaptation utility pickles
with open(indir1 + "adapt_count_all_REF_2000", 'rb') as f:
     days_REF_2000 = cp.load(f)
f.close()

with open(indir1 + "adapt_count_all_REF", 'rb') as f:
     days_REF = cp.load(f)
f.close()

with open(indir1 + "adapt_count_all_P45", 'rb') as f:
     days_P45 = cp.load(f)
f.close()

with open(indir1 + "adapt_count_all_P37", 'rb') as f:
     days_P37 = cp.load(f)
f.close()

# Load other important data
# For REF 2000
# Define list of years
years_2000 = list(range(1981, 2011))

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

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Calculate pc utility by grid cell for 2000
for r, c in rc:
    days_REF_2000[(r, c)] = sum(days_REF_2000[(r, c)].values()) / 30


# Calculate pc utility by grid cell for 2050
days_REF_2050 = {}
days_P45_2050 = {}
days_P37_2050 = {}

for r, c in rc:
    days_REF_2050[(r, c)] = []
    days_P45_2050[(r, c)] = []
    days_P37_2050[(r, c)] = []
    for Yr in year_set[0]:
        for IC in IC_list:
            days_REF_2050[(r, c)].append(days_REF[Yr][IC][(r, c)])
            days_P45_2050[(r, c)].append(days_P45[Yr][IC][(r, c)])
            days_P37_2050[(r, c)].append(days_P37[Yr][IC][(r, c)])
            
    days_REF_2050[(r, c)] = np.mean(days_REF_2050[(r, c)])
    days_P45_2050[(r, c)] = np.mean(days_P45_2050[(r, c)])
    days_P37_2050[(r, c)] = np.mean(days_P37_2050[(r, c)])        
    
# Calculate pc utility by grid cell for 2100
days_REF_2100 = {}
days_P45_2100 = {}
days_P37_2100 = {}

for r, c in rc:
    days_REF_2100[(r, c)] = []
    days_P45_2100[(r, c)] = []
    days_P37_2100[(r, c)] = []
    for Yr in year_set[1]:
        for IC in IC_list:
            days_REF_2100[(r, c)].append(days_REF[Yr][IC][(r, c)])
            days_P45_2100[(r, c)].append(days_P45[Yr][IC][(r, c)])
            days_P37_2100[(r, c)].append(days_P37[Yr][IC][(r, c)])
            
    days_REF_2100[(r, c)] = np.mean(days_REF_2100[(r, c)])
    days_P45_2100[(r, c)] = np.mean(days_P45_2100[(r, c)])
    days_P37_2100[(r, c)] = np.mean(days_P37_2100[(r, c)])    

# Create dataframe of pc utility POL_YEAR vs REF_2000
# Import grid  gdf
grid = gpd.read_file(indir3 + "BenMapGridPoints.csv")
grid = grid.drop('geometry', axis = 1)
grid['geom'] = grid['geom'].apply(wkt.loads)
grid = grid.set_geometry('geom')
grid = grid.set_crs(epsg=4326)
grid['ROW'] = grid['ROW'].astype('int')
grid['COL'] = grid['COL'].astype('int')
grid = grid.drop([1, 219, 220])


# Create columns for each scenario
grid['P45 2050'] = 0
grid['P37 2050'] = 0
grid['P45 2100'] = 0
grid['P37 2100'] = 0

# Fill in dataframe
for r, c in rc:
    grid.loc[(grid['ROW'] == r) & (grid['COL'] == c), 'P45 2050'] = days_P45_2050[(r, c)] - days_REF_2050[(r, c)]
    grid.loc[(grid['ROW'] == r) & (grid['COL'] == c), 'P37 2050'] = days_P37_2050[(r, c)] - days_REF_2050[(r, c)]
    grid.loc[(grid['ROW'] == r) & (grid['COL'] == c), 'P45 2100'] = days_P45_2100[(r, c)] - days_REF_2100[(r, c)]
    grid.loc[(grid['ROW'] == r) & (grid['COL'] == c), 'P37 2100'] = days_P37_2100[(r, c)] - days_REF_2100[(r, c)]
    
    
# Save outputs
with open(outdir + 'adap_days_vs_REFYEAR', 'wb') as f:
      cp.dump(grid,f)
f.close()
    
    
    
    