# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 09:54:11 2022

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
from shapely.geometry import Point

# Load Cell List
grid = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\CAMGridIsaiahBenMap.csv")
grid = grid[grid['In_US'] == True]
grid.reset_index(inplace = True)
grid = grid[['Lon', 'Lat', 'ROW', 'COL']]
grid['geom'] = 0

# Create Shapely points
for i in range(len(grid['Lon'])):
    grid['geom'][i] = Point(grid['Lon'][i].item(),grid['Lat'][i].item())
    
# Drop origional columns
# grid = grid.drop(['Lon', 'Lat'], axis = 1)

# Convert to gdf
grid = gpd.GeoDataFrame(grid) 
grid = grid.set_geometry('geom')
grid.plot()

# Output to pickle and csv
grid.to_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMapGridPoints.csv", index = False)

with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMapGridPoints", 'wb') as f:
     cp.dump(grid,f)
f.close()

