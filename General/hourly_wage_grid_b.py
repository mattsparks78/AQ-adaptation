# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 15:37:17 2022

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

# Import county wage csv file
cwage = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\county_wage.csv", dtype = 'object')

# Import county to grid conversion file
conv = pd.read_pickle(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\county_population_per_cell')

# Import counties with correct names to add data to
cdata = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\county_pop.csv", dtype = 'object', encoding='latin-1')

# connect wage data to cdata
# Create column for wage data in cdata
cdata['hourly_wage'] = 0

# Add wage data based on matching FIPS
for fips in cdata['FIPS']:
    if fips in list(cwage['fips']):
        w = cwage.loc[cwage['fips'] == fips, 'avg_hourly_pay'].item()
        cdata.loc[cdata['FIPS'] == fips, 'hourly_wage'] = w
    else:
        cdata.loc[cdata['FIPS'] == fips, 'hourly_wage'] = 16.30 # national mean wage. Few small counties not included in list

# load AQgrid
with open(r'C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AQGrid\AQGridGeoDF', 'rb') as f:
     AQgrid = cp.load(f)
f.close()
AQgrid['wage'] = 0

# convert county data to grid data
for r, c in zip(list(AQgrid['ROW']),list(AQgrid['COL'])):  
    s = conv[conv['ROW'] == r]
    s = s[s['COL'] == c]
    h = []
    for fips in list(s['FIPS']):
        prop = s.loc[s['FIPS'] == fips, 'pop_share'].item()
        cw = cdata.loc[cdata['FIPS'] == fips, 'hourly_wage'].item()
        h.append(prop * float(cw))
    AQgrid.loc[(AQgrid['ROW'] == r)&(AQgrid['COL'] == c), 'wage'] = sum(h)

# AQgrid.plot('wage')

# Export AQgrid with wage
AQgrid.to_file(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\wage_per_grid.shp')        
AQgrid.to_csv(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\wage_per_grid.csv')


