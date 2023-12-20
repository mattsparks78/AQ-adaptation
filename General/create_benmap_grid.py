# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 09:30:20 2022

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
from shapely import geometry

# Load Cell List
grid = pd.read_csv(r".\Data\CAMGridIsaiahBenMap.csv")
grid = grid[grid['In_US'] == True]
grid.reset_index(inplace = True)
grid = grid[['Lon', 'Lat', 'ROW', 'COL']]

# Calculate cell width/length divided by 2
dLon = 2.5/2
dLat = 1.8947/2

# Lists of Lons and Lats at center of cell
lons = list(grid['Lon'])
lats = list(grid['Lat'])

# Calculate points (-lon, - lat)
# NW point (-lon, + lat)
# NE point (+lon, +lat)
# SE point (+lon, - lat)
SW = []
NW= []
NE = []
SE = []
for i in range(len(lons)):
    a = ((lons[i] - dLon),(lats[i] - dLat))
    SW.append(a)
    b = ((lons[i] - dLon),(lats[i] + dLat))
    NW.append(b)
    c = ((lons[i] + dLon),(lats[i] + dLat))
    NE.append(c)
    d = ((lons[i] + dLon),(lats[i] - dLat))
    SE.append(d)

# Add to dataframe
grid['geom'] = 0

for i in range(len(grid['geom'])):
    # xpt = [SW[i][0], NW[i][0], NE[i][0], SE[i][0]]
    # ypt = [SW[i][1], NW[i][1], NE[i][1], SE[i][1]]
    grid['geom'][i] = geometry.Polygon([SW[i], NW[i], NE[i], SE[i]])
    
# Trim and export
grid = grid[['ROW', 'COL', 'geom']]
grid = gpd.GeoDataFrame(grid)    

grid.to_csv(r".\Data\BenMapGrid.csv", index = False)

with open(r".\Data\BenMapGrid", 'wb') as f:
     cp.dump(grid, f)
     f.close()
