# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:58:20 2022

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

# Import DMA JSON

with open(r".\Data\AQGrid\dma_full.json", 'rb') as f:
    dma_in = json.load(f)
    
f.close()

DMAs = gpd.GeoDataFrame.from_features(dma_in['features'])

# Convert GeometryCollection to MultiPolygon
for i in range(len(DMAs)):
    if DMAs['geometry'][i].type == 'GeometryCollection':
        DMAs['geometry'][i] = unary_union(list((DMAs['geometry'][i])))
    else:
        pass

DMAs.crs = "EPSG:4326"

# Import AQgrid
with open(r'.\Data\AQGrid\AQGridGeoDF', 'rb') as f:
     AQgrid = cp.load(f)
f.close()

AQgrid.crs = "EPSG:4326"

# Reproject to Albers to calculate more accurate areas
DMAs = DMAs.to_crs("EPSG:5070")
AQgrid = AQgrid.to_crs("EPSG:5070")

# Compute intersections between DMA grid and AQgrid
# Calculate areas for grids
DMAs['area_DMA'] = DMAs.area
AQgrid['area_grid'] = AQgrid.area

# Find intersections
joined = gpd.overlay(AQgrid, DMAs, how = 'intersection')

# Calculate Joined area
joined['area_joined'] = joined.area

# Caclulcate percentage of DMA area in each grid cell
joined['percentage_area'] = (joined['area_joined']/joined['area_DMA'])

# Output the important information
harmonized = joined[['ROW', 'COL', "NAME", 'DMA', 'percentage_area']].copy()

# Pickle harmonized
harmonized.to_pickle(r'.\Data\dma_per_cell')

# Test
test = pd.read_pickle(r'.\Data\dma_per_cell')
