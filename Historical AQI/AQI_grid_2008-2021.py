# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 16:10:15 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
import scipy.stats
from statsmodels.distributions.empirical_distribution import ECDF

# Import grid gdf
with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\BenMapGrid", 'rb') as f:
     AQgrid = cp.load(f)
f.close()

AQgrid = AQgrid.drop('geom', axis = 1)

# Define list of grid cells
rc = list(zip(AQgrid['ROW'],AQgrid['COL']))

# Set Directories
indir1 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\EPA Standard Monitor Data\\'
indir2 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMAP Neighbor Files\\'
outdir = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical AQI\\'

# Create dict to hold PM and O3 data.
o3_dict = {}
pm_dict = {}

# Initialize Years
for Yr in range(2008, 2022):
    o3_dict[Yr] = {}
    pm_dict[Yr] = {}


# Convert PM and O3 sensor data to grid data
for Yr in range(2008,2022):
    # Create filenames
    fm_pm = 'EPA Standard Monitors_PM2.5_' + str(Yr) +'.csv'
    fm_o = 'EPA Standard Monitors_Ozone_' + str(Yr) + '.csv'
    fn_pm = 'PM_' + str(Yr) + '_neighborfile.csv'
    fn_o3 = 'O3_' + str(Yr) + '_neighborfile.csv'
    
    # Load dataframes
    mon_pm = pd.read_csv(indir1 + fm_pm)
    mon_o3 = pd.read_csv(indir1 + fm_o)
    nei_pm = pd.read_csv(indir2 + fn_pm)
    nei_o3 = pd.read_csv(indir2 + fn_o3)
    
    # Convert to dict
    mon_o3 = mon_o3[['Monitor Name', 'Values']]
    mon_o3_dict = dict(mon_o3.values)
    mon_pm = mon_pm[['Monitor Name', 'Values']]
    mon_pm_dict = dict(mon_pm.values)
    
    
    for poll in ['O3', 'PM25']:
        if poll == 'O3':
            for r, c in rc:
                # Isolate sensor weights for current grid cell
                n = nei_o3[(nei_o3['Row'] == r) & (nei_o3['Col'] == c)]
                # Reset index to make easier to work with
                n = n.reset_index()
                # Create storage list for each relevant sensor
                o = [None] * len(n)
                # Create new list for output data
                o3_out = []
                # Create key to index output 
                key = '(' + str(r) + ', ' + str(c) + ")"
                # Iterate over each sensor 
                for i in range(len(n)):
                    # Get current sensor name
                    mon = n['MonitorName'][i]
                    # Get raw data from sensor
                    o[i] = mon_o3_dict[mon]
                    # Split from string to list
                    o[i] = o[i].split(',')
                    # Change missing entries to nan
                    o[i] = [np.nan if x == '.' else x for x in o[i]]
                    # Change measured values from str to int
                    o[i] = [int(x) if type(x) == str else x for x in o[i]]
                    # convert list to np array
                    o[i] = np.asarray(o[i])
                # Iterate over each day within o
                for j in range(len(o[0])):
                    # Create list for sensor data/distances
                    conc = []
                    dist = []
                    for k in range(len(o)):
                    # Determine which data we have for day j
                        if o[k][j] >= 0:
                            conc.append(o[k][j])
                            dist.append(n['Distance'][k])
                        else:
                            pass
                    # Initialize storage for 1/d, weights, and weighted concentrations
                    inv_dist = []
                    weights = []
                    weighted_con = []
                    #convert distances to weights
                    #
                    for l in range(len(dist)):
                        inv_dist.append((1/(dist[l])))
                    sum_inv = sum(inv_dist)
                    for l in range(len(dist)):
                        weights.append((inv_dist[l])/sum_inv)
                    for i in range(len(conc)):
                        weighted_con.append(round(conc[i] * weights[i]))
                    o3_out.append(sum(weighted_con))        
                o3_dict[Yr][key] = o3_out
        else:
            for r, c in rc:
                # Isolate sensor weights for current grid cell
                n = nei_pm[(nei_pm['Row'] == r) & (nei_pm['Col'] == c)]
                # Ensure that sensors in neighbor file correspond with measurements
                for mon in n['MonitorName']:
                    if mon in list(mon_pm_dict.keys()):
                        pass
                    else:
                        n = n.drop(n[n.MonitorName == mon].index)
                if len(n) == 0:
                    pm_dict[Yr][key] = list(0 for i in range(len(o3_out)))
                else:
                    # Reset index to make easier to work with
                    n = n.reset_index()
                    # Create storage list for each relevant sensor
                    o = [None] * len(n)
                    # Create new list for output data
                    pm_out = []
                    # Create key to index output 
                    key = '(' + str(r) + ', ' + str(c) + ")"
                    # Iterate over each sensor 
                    for i in range(len(n)):
                        # Get current sensor name
                        mon = n['MonitorName'][i]
                        # Get raw data from sensor
                        o[i] = mon_pm_dict[mon]
                        # Split from string to list
                        o[i] = o[i].split(',')
                        # Change missing entries to nan
                        o[i] = [np.nan if x == '.' else x for x in o[i]]
                        # Change measured values from str to int
                        o[i] = [float(x) if type(x) == str else x for x in o[i]]
                        # convert list to np array
                        o[i] = np.asarray(o[i])
                    # Iterate over each day within o
                    for j in range(len(o[0])):
                        # Create list for sensor data/distances
                        conc = []
                        dist = []
                        for k in range(len(o)):
                        # Determine which data we have for day j
                            if o[k][j] >= 0:
                                conc.append(o[k][j])
                                dist.append(n['Distance'][k])
                            else:
                                pass
                        # Initialize storage for 1/d, weights, and weighted concentrations
                        inv_dist = []
                        weights = []
                        weighted_con = []
                        # convert distances to weights
                        # Calculate 1/distance
                        for l in range(len(dist)):
                            inv_dist.append((1/(dist[l])))
                        # Calculate sum of 1/distances for denominator
                        sum_inv = sum(inv_dist)
                        # Calculate weights for each sensor
                        for l in range(len(dist)):
                            weights.append((inv_dist[l])/sum_inv)
                        # Calculate contribution from each sensor
                        for i in range(len(conc)):
                            weighted_con.append(round(conc[i] * weights[i]))
                        pm_out.append(sum(weighted_con))        
                    # Output daily measurements for grid cell
                    pm_dict[Yr][key] = pm_out

# Output pm_dict, o3_dict pickles
with open(outdir + 'hist_pm_2008-2021', 'wb') as f:
      cp.dump(pm_dict,f)
f.close()

with open(outdir + 'hist_o3_2008-2021', 'wb') as f:
      cp.dump(o3_dict,f)
f.close()

# Calculate AQI based on pm and o3 measurements
# Create dicts for O3 and PM AQI
pm_aqi = copy.deepcopy(pm_dict)
o3_aqi = copy.deepcopy(o3_dict)

# Iterate over years, grid cells, days
for Yr in range(2008,2022):
    for cell in pm_aqi[Yr]:
        # Iterate over items in array
        for i in range(len(pm_aqi[Yr][cell])):
            con = pm_dict[Yr][cell][i]
            if con <= 12:
                aqi = round(4.17*con)
            elif 12 < con <= 35.4:
                aqi = round(2.1*(con - 12) + 51)
            elif 35.4 < con <= 55.4:
                aqi = round(2.46*(con - 35.4) + 101)
            elif 55.4 < con <= 150.4:
                aqi = round(0.52*(con-55.4) + 151)
            elif 150.4 < con <= 250.4:
                aqi = round(0.99*(con-150.4) + 201)
            else:
                aqi = round(0.8*(con-250.4) + 301)
            pm_aqi[Yr][cell][i] = aqi
    
    for cell in o3_aqi[Yr]:
        # Iterate over items in array
        for i in range(len(o3_aqi[Yr][cell])):
            con = o3_dict[Yr][cell][i]
            if con <= 54:
                aqi = round(0.93*con)
            elif 54 < con <= 70:
                aqi = round(3.06*(con-54) + 51)
            elif 70 < con <= 85:
                aqi = round(3.27*(con - 70)+ 101)
            elif 85 < con <= 105:
                aqi = round(2.58*(con - 85) + 151)
            elif 105 < con <= 200:
                aqi = round(1.05*(con - 105) + 201)
            else:
                aqi = 301
            o3_aqi[Yr][cell][i] = aqi

# Create dict to hold results
aqi_dict = copy.deepcopy(pm_aqi)            
for Yr in range(2008,2022):
    for cell in aqi_dict[Yr]:
        for i in range(len(pm_aqi[Yr][cell])):
            pm = pm_aqi[Yr][cell][i]
            o3 = o3_aqi[Yr][cell][i]
            val = max(pm,o3)
            aqi_dict[Yr][cell][i] = val

# Output aqi_dict as pickle
with open(outdir + 'hist_aqi_2008-2021', 'wb') as f:
      cp.dump(aqi_dict,f)
f.close()

with open(outdir + 'hist_aqi_2008-2021_PM', 'wb') as f:
      cp.dump(pm_aqi,f)
f.close()

with open(outdir + 'hist_aqi_2008-2021_O3', 'wb') as f:
      cp.dump(o3_aqi,f)
f.close()
