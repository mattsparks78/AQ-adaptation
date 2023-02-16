# -*- coding: utf-8 -*-
"""
Created on Thu May  5 12:50:46 2022

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
from geopandas.tools import overlay

# Import grid gdf
with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\BenMapGrid", 'rb') as f:
     AQgrid = cp.load(f)
f.close()

# Import county shp file
counties_full = gpd.read_file(r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AQGrid\co99_d00.shp")

# Delete HI, AK, PR, Guam
counties_full.drop(counties_full.loc[counties_full['STATE']== '72'].index, inplace=True)
counties_full.drop(counties_full.loc[counties_full['STATE']== '02'].index, inplace=True)
counties_full.drop(counties_full.loc[counties_full['STATE']== '15'].index, inplace=True)

counties = counties_full.copy(deep=True)

# counties gdf has multiple entries for each county. Need to aggregate together
# Create new gdf to hold aggregated counties data
c = gpd.GeoDataFrame()

# Create list of states
states_list = sorted(list(set(counties['STATE'])))

# Aggregate county geometries
for state in states_list:
    temp_state = gpd.GeoDataFrame()
    temp_state = (counties[counties['STATE'] == state])
    temp_merged = temp_state.dissolve(by = 'COUNTY')
    c = c.append(temp_merged)



# tests to see if combined correctly
# c.boundary.plot()
# counties.boundary.plot()

# sum(c['AREA'])
# sum(counties['AREA'])

# Determine proportion of county within each grid cell
# Set CRS for each gdf
c = c.set_crs("EPSG:4326")
AQgrid = AQgrid.set_geometry('geom')
AQgrid = AQgrid.set_crs("EPSG:4326")

# Pull out county number data from index
c['COUNTY_CODE'] = c.index

# Combine state code and county code
c['FIPS'] = c['STATE'] + c['COUNTY_CODE']

# Pull out important information from c
ci = c[['NAME', 'FIPS', 'geometry']]

# Compute intersections between county and AQgrid
# Reproject to Albers to calculate area more accurately
ci = ci.to_crs("EPSG:5070")
AQgrid = AQgrid.to_crs("EPSG:5070")

# Calculate areas for grids
ci['area_county'] = ci.area
AQgrid['area_grid'] = AQgrid.area

# Find intersections
joined = gpd.overlay(AQgrid, ci, how = 'intersection')

# Calculate Joined area
joined['area_joined'] = joined.area

# Caclulcate percentage of DMA area in each grid cell
joined['percentage_area'] = (joined['area_joined']/joined['area_county'])

# Output the important information
harmonized = joined[['ROW', 'COL', "NAME", 'FIPS', 'percentage_area']].copy()

# Pickle harmonized
# harmonized.to_csv(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\county_per_cell_b.csv', index = False)

# Test
# test = pd.read_pickle(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\county_per_cell_b')
