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
with open(r".\Data\BenMap Grid\BenMapGrid", 'rb') as f:
     AQgrid = cp.load(f)
f.close()

AQgrid = AQgrid.drop([1, 219, 220], axis = 0)
AQgrid = AQgrid.drop('geom', axis = 1)

# Define list of grid cells
rc = list(zip(AQgrid['ROW'],AQgrid['COL']))

# Set Directories
indir1 = r'.\Data\EPA Standard Monitor Data\\'
indir2 = r'.\Data\BenMAP Neighbor Files\\'
indir3 = r'.\Data\General\\'
outdir = r'.\Data\Historical AQI\\'

# Create dict to hold PM and O3 data.
o3_dict = {}
pm_dict = {}

# Initialize Years
for Yr in range(2013, 2022):
    o3_dict[Yr] = {}
    pm_dict[Yr] = {}


# Convert PM and O3 sensor data to grid data
for Yr in range(2013,2022):
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
                    pm_dict[Yr][(r, c)] = pm_out

# Calculate mortalities based on pm_dict
# Load necessary data
# Import hourly wage data
wage = pd.read_csv(indir3 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Import population per grid cell
pop = pd.read_csv(indir3 + "pop_lepeule.csv")

# Import baseline mortality per grid cell
bm = pd.read_csv(indir3 + 'baseline_mortality_b.csv')

# Set RR 1.06 for Krewski, 1.14 for Lepeule
RR = 1.14
AF = (RR - 1)/ RR

# Set VSL
VSL = 7400000

# Calc mort and valuation
test = copy.deepcopy(pm_dict)
test[2013][(9, 16)] = 0
test[2013][(12, 24)] = 0
test[2013][(13, 25)] = 0
mean_pm_cells = {}

# Drop 0 values from pm_dict
for Yr in range(2013,2022):
    mean_pm_cells[Yr] = {}
    for r, c in rc:
         vals = test[Yr][(r, c)]
         if np.sum(vals) != 0:
             vals = list(filter(lambda num: num != 0, vals))
         else:
             vals = 0
         mean_pm_cells[Yr][(r, c)] = np.mean(vals)

mean_pm_cells[2013][(9, 16)] = mean_pm_cells[2014][(9, 16)]
mean_pm_cells[2013][(12, 24)] = mean_pm_cells[2014][(12, 24)]
mean_pm_cells[2013][(13, 25)] = mean_pm_cells[2014][(13, 25)]

# Calculate full burden mortalities based on pm concentrations
pm_mort_count = {}
pm_mort_val = {}
for Yr in range(2013, 2022):
    pm_mort_count[Yr] = {}
    pm_mort_val[Yr] = {}
    for r, c in rc:
        popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item()
        bmr = bm.loc[(bm['ROW'] == r) & (bm['COL'] == c), 'BaseMortRate'].item()
        dPM = mean_pm_cells[Yr][(r, c)]
        pm_mort_count[Yr][(r, c)] = bmr * AF * dPM * 0.1 * popu
        pm_mort_val[Yr][(r, c)] = pm_mort_count[Yr][(r, c)] * VSL

# Calculate national valuations
pm_nat_val = {}
for Yr in range(2013, 2022):
    pm_nat_val[Yr] = []
    for r, c in rc:
        pm_nat_val[Yr].append(pm_mort_val[Yr][(r, c)])
       
    pm_nat_val[Yr] = np.sum(pm_nat_val[Yr])    

# Output pickles
# with open(outdir + 'morbidity_count_2013-2021', 'wb') as f:
#      cp.dump(pm_mort_count,f)
# f.close()

# with open(outdir + 'morbidity_val_by_cell_2013-2021', 'wb') as f:
#      cp.dump(pm_mort_val,f)
# f.close()

# with open(outdir + 'morbidity_val_nation_2013-2021', 'wb') as f:
#      cp.dump(pm_nat_val,f)
# f.close()


