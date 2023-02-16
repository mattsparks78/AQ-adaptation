# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:57:36 2022

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

# Import county AQI data from EPA
aqi_15 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2015.csv")
aqi_16 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2016.csv")
aqi_17 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2017.csv")
aqi_18 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2018.csv")
aqi_19 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2019.csv")
aqi_20 = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA County AQI\daily_aqi_by_county_2020.csv")

# Import county to grid cell conversion df
ctg = pd.read_pickle(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\county_population_per_cell")

# Import base AQ grid to add columns to for each year
gbad = pd.read_pickle(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\CAM_grid_pop')
gbad = gbad.drop(['geom', 'pop'], axis = 1)

gbad15 = gbad.copy()
gbad16 = gbad.copy()
gbad17 = gbad.copy()
gbad18 = gbad.copy()
gbad19 = gbad.copy()
gbad20 = gbad.copy()
# Add FIPS column to AQI dataframes
aqi_15 = aqi_15.astype({'State Code': str, 'County Code': str}, copy = False)
states15 = list(aqi_15['State Code'])
counties15 = list(aqi_15['County Code'])
fips15 = counties15.copy()

for i in range(len(states15)):
    if len(states15[i]) == 1:
        states15[i] = '0' + states15[i]
    else:
        pass
    if len(counties15[i]) == 1:
        counties15[i] = '00' + counties15[i]
    elif len(counties15[i]) == 2:
        counties15[i] = '0' + counties15[i]
    else:
        pass
    fips15[i] = states15[i] + counties15[i]
    
aqi_15.insert(10, 'fips', pd.Series(fips15), allow_duplicates= True)
aqi_15 = aqi_15[['fips', 'AQI', 'Date']]

# 2016
aqi_16 = aqi_16.astype({'State Code': str, 'County Code': str}, copy = False)
states16 = list(aqi_16['State Code'])
counties16 = list(aqi_16['County Code'])
fips16 = counties16.copy()

for i in range(len(states16)):
    if len(states16[i]) == 1:
        states16[i] = '0' + states16[i]
    else:
        pass
    if len(counties16[i]) == 1:
        counties16[i] = '00' + counties16[i]
    elif len(counties16[i]) == 2:
        counties16[i] = '0' + counties16[i]
    else:
        pass
    fips16[i] = states16[i] + counties16[i]
    
aqi_16.insert(10, 'fips', pd.Series(fips16), allow_duplicates= True)
aqi_16 = aqi_16[['fips', 'AQI', 'Date']]

#2017
aqi_17 = aqi_17.astype({'State Code': str, 'County Code': str}, copy = False)
states17 = list(aqi_17['State Code'])
counties17 = list(aqi_17['County Code'])
fips17 = counties17.copy()

for i in range(len(states17)):
    if len(states17[i]) == 1:
        states17[i] = '0' + states17[i]
    else:
        pass
    if len(counties17[i]) == 1:
        counties17[i] = '00' + counties17[i]
    elif len(counties17[i]) == 2:
        counties17[i] = '0' + counties17[i]
    else:
        pass
    fips17[i] = states17[i] + counties17[i]
    
aqi_17.insert(10, 'fips', pd.Series(fips17), allow_duplicates= True)
aqi_17 = aqi_17[['fips', 'AQI', 'Date']]

#2018
aqi_18 = aqi_18.astype({'State Code': str, 'County Code': str}, copy = False)
states18 = list(aqi_18['State Code'])
counties18 = list(aqi_18['County Code'])
fips18 = counties18.copy()

for i in range(len(states18)):
    if len(states18[i]) == 1:
        states18[i] = '0' + states18[i]
    else:
        pass
    if len(counties18[i]) == 1:
        counties18[i] = '00' + counties18[i]
    elif len(counties18[i]) == 2:
        counties18[i] = '0' + counties18[i]
    else:
        pass
    fips18[i] = states18[i] + counties18[i]
    
aqi_18.insert(10, 'fips', pd.Series(fips18), allow_duplicates= True)
aqi_18 = aqi_18[['fips', 'AQI', 'Date']]

#2019
aqi_19 = aqi_19.astype({'State Code': str, 'County Code': str}, copy = False)
states19 = list(aqi_19['State Code'])
counties19 = list(aqi_19['County Code'])
fips19 = counties19.copy()

for i in range(len(states19)):
    if len(states19[i]) == 1:
        states19[i] = '0' + states19[i]
    else:
        pass
    if len(counties19[i]) == 1:
        counties19[i] = '00' + counties19[i]
    elif len(counties19[i]) == 2:
        counties19[i] = '0' + counties19[i]
    else:
        pass
    fips19[i] = states19[i] + counties19[i]
    
aqi_19.insert(10, 'fips', pd.Series(fips19), allow_duplicates= True)
aqi_19 = aqi_19[['fips', 'AQI', 'Date']]

#2020
aqi_20 = aqi_20.astype({'State Code': str, 'County Code': str}, copy = False)
states20 = list(aqi_20['State Code'])
counties20 = list(aqi_20['County Code'])
fips20 = counties20.copy()

for i in range(len(states20)):
    if len(states20[i]) == 1:
        states20[i] = '0' + states20[i]
    else:
        pass
    if len(counties20[i]) == 1:
        counties20[i] = '00' + counties20[i]
    elif len(counties20[i]) == 2:
        counties20[i] = '0' + counties20[i]
    else:
        pass
    fips20[i] = states20[i] + counties20[i]
    
aqi_20.insert(10, 'fips', pd.Series(fips20), allow_duplicates= True)
aqi_20 = aqi_20[['fips', 'AQI', 'Date']]

# Add columns for each time period of interest
gbad15 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])
gbad16 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])
gbad17 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])
gbad18 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])
gbad19 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])
gbad20 = gbad15.reindex(columns = ['ROW', 'COL', 'JAN', 'FEB', 'MAR','APR','MAY','JUN','JUL', 'AUG','SEP','OCT','NOV','DEC'])    

# Define rc
rc = list(zip(list(gbad['ROW']),list(gbad['COL'])))  

# Identify bad air days for each grid cell
for r, c in rc:
    counties = ctg[(ctg['ROW'] == r) & (ctg['COL'] == c)]
    for fips in counties['FIPS']:
        if fips in fips15:
            m = aqi_15.loc[aqi_15['fips'] == str(fips)]
            bad = m.loc[m['AQI'] >= 100]
            
        else:
            pass