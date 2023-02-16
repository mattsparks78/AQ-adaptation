# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 12:56:06 2022

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

indir1 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\'
indir2 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\US Boundary\\'
outdir = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\'

# Load BenMap Grid
with open(indir1 + 'BenMapGrid', 'rb') as f:
     BMgrid = cp.load(f)
f.close()

BMgrid.drop([1, 219, 220], axis = 0, inplace = True)
BMgrid = BMgrid.reset_index(drop = True)
BMgrid = BMgrid.set_geometry('geom')
BMgrid = BMgrid.set_crs('epsg:4326')

# Create columns for proportion of cell in each census region
grid_area = copy.deepcopy(BMgrid)

grid_area['Northeast'] = 0
grid_area['Midwest'] = 0
grid_area['South'] = 0
grid_area['West'] = 0

grid_prop = copy.deepcopy(grid_area)

# Create shapefiles for each census region
usbound = gpd.read_file(indir2 + 'tl_2012_us_state.shp')

usbound = usbound.to_crs('epsg:4326')


northeast = usbound.loc[usbound['REGION'] == '1']
midwest = usbound.loc[usbound['REGION'] == '2']
south = usbound.loc[usbound['REGION'] == '3']
west = usbound.loc[usbound['REGION'] == '4']
# Drop HI and AK
west.drop([0, 47], axis = 0, inplace = True)

# Calculate intersection of grid cells and census regions
merge_ne = gpd.overlay(grid_area, northeast, how = 'intersection')
merge_mw = gpd.overlay(grid_area, midwest, how = 'intersection')
merge_s = gpd.overlay(grid_area, south, how = 'intersection')
merge_w = gpd.overlay(grid_area, west, how = 'intersection')

merge_ne.reset_index(drop = True, inplace = True)
merge_mw.reset_index(drop = True, inplace = True)
merge_s.reset_index(drop = True, inplace = True)
merge_w.reset_index(drop = True, inplace = True)

for i in range(len(merge_ne)):
    s = merge_ne.iloc[i]
    s = pd.DataFrame(s)
    s = s.transpose()
    r = s['ROW'].item()
    c = s['COL'].item()
    area = s['Shape_Area'].item()
    
    grid_area.loc[(grid_area['ROW'] == r) & (grid_area['COL'] == c), 'Northeast'] += area

for i in range(len(merge_mw)):
    s = merge_mw.iloc[i]
    s = pd.DataFrame(s)
    s = s.transpose()
    r = s['ROW'].item()
    c = s['COL'].item()
    area = s['Shape_Area'].item()
    
    grid_area.loc[(grid_area['ROW'] == r) & (grid_area['COL'] == c), 'Midwest'] += area

for i in range(len(merge_s)):
    s = merge_s.iloc[i]
    s = pd.DataFrame(s)
    s = s.transpose()
    r = s['ROW'].item()
    c = s['COL'].item()
    area = s['Shape_Area'].item()
    
    grid_area.loc[(grid_area['ROW'] == r) & (grid_area['COL'] == c), 'South'] += area

for i in range(len(merge_w)):
    s = merge_w.iloc[i]
    s = pd.DataFrame(s)
    s = s.transpose()
    r = s['ROW'].item()
    c = s['COL'].item()
    area = s['Shape_Area'].item()
    
    grid_area.loc[(grid_area['ROW'] == r) & (grid_area['COL'] == c), 'West'] += area


# Calculate proportion of grid cell in each census region
for i in range(len(grid_prop)):
    t = grid_area.iloc[i]
    t = pd.DataFrame(t)
    t = t.transpose()
    total_area = t['Northeast'].item() + t['Midwest'].item() + t['South'].item() + t['West'].item()
    r = t['ROW'].item()
    c = t['COL'].item()
    
    grid_prop.loc[(grid_prop['ROW'] == r) & (grid_prop['COL'] == c), 'Northeast'] = t['Northeast'].item() / total_area
    grid_prop.loc[(grid_prop['ROW'] == r) & (grid_prop['COL'] == c), 'Midwest'] = t['Midwest'].item() / total_area
    grid_prop.loc[(grid_prop['ROW'] == r) & (grid_prop['COL'] == c), 'South'] = t['South'].item() / total_area
    grid_prop.loc[(grid_prop['ROW'] == r) & (grid_prop['COL'] == c), 'West'] = t['West'].item() / total_area
    
# Output files
# with open(outdir + 'grid_by_census_region', 'wb') as f:
#      cp.dump(grid_prop,f)
# f.close()

# grid_prop.to_csv(outdir + 'grid_by_census_region.csv', index = False)

west_merge = west['geometry'].unary_union
ne_merge = northeast['geometry'].unary_union
s_merge = south['geometry'].unary_union
mw_merge = midwest['geometry'].unary_union

list_regs = ['Northeast', 'Midwest', 'South', 'West']

list_geoms = [ne_merge, mw_merge, s_merge, west_merge]


region_geoms = gpd.GeoDataFrame(data = list_regs, geometry = list_geoms)
region_geoms = region_geoms.rename(columns = {0: 'regions'})


region_geoms.to_file(outdir + "US_census_regions.shp")
