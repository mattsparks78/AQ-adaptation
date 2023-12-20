# -*- coding: utf-8 -*-
"""
Created on Fri May 13 10:50:05 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
import json
from geopandas.tools import overlay

# Load AQgrid
with open(r'.\Data\BenMap Grid\BenMapGrid', 'rb') as f:
     AQgrid = cp.load(f)
f.close()

AQgrid.crs = "EPSG:4326"

# Load Counties to AQgrid
counties_grid = pd.read_csv(r'.\Data\BenMap Grid\county_per_cell_b.csv', dtype = 'object')
counties_grid = counties_grid.astype({'ROW':int, 'COL':int, 'FIPS': str, 'percentage_area': float})
# Add population data to counties #ADD BEDFORD CITY TO CSV FILE
counties_pop = pd.read_csv(r".\Data\General\county_pop.csv", dtype = 'object', encoding='latin-1')
# Convert population data back to int. Brought all in as obj to keep leading 0 in FIPS
counties_pop = counties_pop.astype({'POPESTIMATE2019': int})


# Create list of row,column pairs in grid cell
rc = list(zip(list(AQgrid['ROW']),list(AQgrid['COL'])))

# Create new df with list of grid cells in one column and a corresponding column for each grid cell
sci_grid = pd.DataFrame(data = None, columns = rc)
sci_grid.insert(0,'base_cell', rc)

# Calculate population share for each county in each cell
counties_grid['pop'] = 0



# Calculate population within each county in grid cell
for i in range(len(counties_grid)):
    fips = counties_grid['FIPS'][i]
    tot_pop = counties_pop.loc[counties_pop['FIPS'] == fips, 'POPESTIMATE2019'].item()
    counties_grid['pop'][i] = round(counties_grid['percentage_area'][i]*tot_pop)

# Create column for population share for each county in each cell
counties_grid['pop_share'] = 0

# Calculate grid cell population share by county
# Create data frame that samples the data for the cell corresponding with r, c
for r, c in zip(list(AQgrid['ROW']),list(AQgrid['COL'])):  
    s = counties_grid[counties_grid['ROW'] == r]
    s = s[s['COL'] == c]
    s['pop_share'] = s['pop']/sum(s['pop'])
    
    for fips_out in list(s['FIPS']):
            counties_grid.loc[(counties_grid['ROW'] == r) & (counties_grid['COL'] == c) & (counties_grid['FIPS'] == fips_out), 'pop_share'] = s.loc[s['FIPS'] == fips_out, 'pop_share'].item()
    
# save counties_grid
counties_grid.to_pickle(r'.\Data\BenMap Grid\county_population_per_cell_b')
counties_grid.to_csv(r'.\Data\BenMap Grid\county_population_per_cell_b.csv', index = False)


