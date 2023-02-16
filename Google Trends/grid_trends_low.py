# -*- coding: utf-8 -*-
"""
Created on Tue May 10 13:47:33 2022

@author: MSS
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
import json
from shapely.ops import unary_union
from geopandas.tools import overlay


# Create time series of google trends data organized by grid cell

# Load AQgrid
with open(r'C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AQGrid\AQGridGeoDF', 'rb') as f:
     AQgrid = cp.load(f)
f.close()

# Set coordinate reference system
AQgrid.crs = "EPSG:4326"

# Load DMA to AQgrid
dma_grid = pd.read_pickle(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\dma_per_cell')

# Load trends by DMA
dma_trends = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\trends_by_dma.csv")

# Create list of DMAs
dma_list = sorted(list(set(dma_grid['DMA'])))

# Aggregate trends data from dma to grid cells
# Add date columns to AQgrid df
for i in range(len(dma_trends)):
    col_name = dma_trends['date'][i]
    AQgrid[col_name] = 0

# Artifically set to 0. Did not originally have data Glendive DMA which is the smallest (~4300 pop) and covers rural Montana
dma_trends['798_max_ratio_lo'] = 0 

# Create loop to run through all grid cells
for r, c in zip(list(AQgrid['ROW']),list(AQgrid['COL'])):  
# Create data frame that samples the data for the cell corresponding with r, c
    s = dma_grid[dma_grid['ROW'] == r]
    s = s[s['COL'] == c]
# Create list of DMA regions included in cell    
    l = list(s['DMA'])
# Create dataframe to temporarily hold trends data
    ar = pd.DataFrame() 
    ar['date'] = dma_trends['date']
    ar['max_ratio_lo'] = 0.0
        
    for dma in l:
        p = s[s['DMA'] == dma]['percentage_area'].item() 
        ar[dma+'_max_ratio_lo'] = (dma_trends[dma+'_max_ratio_lo'])*p
        ar['max_ratio_lo'] += ar[dma+'_max_ratio_lo']
        ar.drop(dma+'_max_ratio_lo', axis = 1, inplace = True)
    
    
    for d in list(ar['date']):
        AQgrid.loc[(AQgrid['ROW'] == r) & (AQgrid['COL'] == c), d] = ar.loc[ar['date'] == d, 'max_ratio_lo'].item()
 
# Output shapefile
AQgrid.to_file(driver = 'ESRI Shapefile', filename = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\trends_cell_lo.shp')

# Test input pickle
test = gpd.read_file(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\trends_cell_lo.shp')

# Output csv
AQgrid.to_csv(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\trends_cell_lo.csv')

# Plot sample data
fig, ax = plt.subplots(1)
ax = AQgrid.plot(column = '2017-01-22', ax = ax, legend = True)
fig.suptitle('Google Trends Data by Grid 2017-01-22')
plt.show()

fig.savefig(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\Trends_lo_2017_01_22', dpi = 300)
