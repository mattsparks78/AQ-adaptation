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

AQgrid = AQgrid.drop('geom', axis = 1)

# Set Directories
indir1 = r'.\Data\Historical Model PM\\'
indir2 = r'.\Data\Historical Model O3\\'
outdir = r'.\Data\Historical Model PM\\'

# Load pickles of pm, o3 data
with open(indir1 + 'pm_conc_base_cell_model_2000', 'rb') as f:
     pm_dict = cp.load(f)
f.close()

with open(indir2 + 'o3_conc_base_cell_model_2000', 'rb') as f:
     o3_dict = cp.load(f)
f.close()

# Define list of grid cells
rc = list(pm_dict.keys())

# Initialize list to hold all data
pm_list = []
o3_list = []

# Add pollutant concentrations to lists
for r, c in rc:
    for Yr in range(1981,2011):
        pm_list.extend(pm_dict[(r, c)][Yr])
        o3_list.extend(o3_dict[(r, c)][Yr])

# Create lists for aqi data
pm_aqi_list = []
o3_aqi_list = []

# Iterate over length of lists
for i in range(len(pm_list)):
    con_pm = pm_list[i]
    con_o3 = o3_list[i]
    
    if con_pm <= 12:
        aqi_pm = (4.17*con_pm)
    elif 12 < con_pm <= 35.4:
        aqi_pm = (2.1*(con_pm - 12) + 51)
    elif 35.4 < con_pm <= 55.4:
        aqi_pm = (2.46*(con_pm - 35.4) + 101)
    elif 55.4 < con_pm <= 150.4:
        aqi_pm = (0.52*(con_pm-55.4) + 151)
    elif 150.4 < con_pm <= 250.4:
        aqi_pm = (0.99*(con_pm-150.4) + 201)
    else:
        aqi_pm = (0.8*(con_pm-250.4) + 301)
    
    pm_aqi_list.append(aqi_pm)
    
    if con_o3 <= 54:
        aqi_o3 = (0.93*con_o3)
    elif 54 < con_o3 <= 70:
        aqi_o3 = (3.06*(con_o3-54) + 51)
    elif 70 < con_o3 <= 85:
        aqi_o3 = (3.27*(con_o3 - 70)+ 101)
    elif 85 < con_o3 <= 105:
        aqi_o3 = (2.58*(con_o3 - 85) + 151)
    elif 105 < con_o3 <= 200:
        aqi_o3 = (1.05*(con_o3 - 105) + 201)
    else:
        aqi_o3 = 301
    
    o3_aqi_list.append(aqi_o3)
    
tot_aqi_list = []

for i in range(len(pm_list)):
    tot_aqi_list.append(max(pm_aqi_list[i], o3_aqi_list[i]))

dist = ECDF(tot_aqi_list)

pctl_100 = dist(100)
pctl_150 = dist(150)

prop_o3_100 = len([i for i in o3_aqi_list if i > 100]) / (len([i for i in o3_aqi_list if i > 100]) + len([i for i in pm_aqi_list if i > 100]))
prop_pm_100 = len([i for i in pm_aqi_list if i > 100]) / (len([i for i in pm_aqi_list if i > 100]) + len([i for i in o3_aqi_list if i > 100]))

prop_o3_150 = len([i for i in o3_aqi_list if i > 150]) / (len([i for i in o3_aqi_list if i > 150]) + len([i for i in pm_aqi_list if i > 150]))
prop_pm_150 = len([i for i in pm_aqi_list if i > 150]) / (len([i for i in pm_aqi_list if i > 150]) + len([i for i in o3_aqi_list if i > 150]))


plt.hist(o3_aqi_list, bins = 50, density = True, histtype = 'step', stacked = True, fill = False)
plt.hist(pm_aqi_list, bins = 50, density = True, histtype = 'step', stacked = True, fill = False)
plt.xlabel('AQI')
plt.ylabel('Probability')
plt.title('Modeled Histogram of AQI Probability by Pollutant')
plt.xlim(0, 150)
plt.legend(['O3', 'PM'])
# plt.savefig(outdir + 'Modeled Histogram by pollutant.png', dpi = 400, bbox_inches = 'tight')

# Plot measured vs modeled pm (with calc_AQI_150_pctl_cell_nat_meas_mod)... what a file name lol
with open(indir1 + 'measured_pm_aqi_2008-2021', 'rb') as f:
     pm_list_meas = cp.load(f)
f.close()

with open(indir2 + 'measured_o3_aqi_2008-2021', 'rb') as f:
     o3_list_meas = cp.load(f)
f.close()


# Plot PM model vs measured
plt.hist(pm_aqi_list, bins = 50, density = True, histtype = 'step', stacked = True, fill = False) # Modeled
plt.hist(pm_list_meas, bins = 50, density = True, histtype = 'step', stacked = True, fill = False) # Measured
plt.xlabel('AQI')
plt.ylabel('Probability')
plt.title('Measured vs Modeled Historical AQI Concentration based on PM')
plt.xlim(0, 150)
plt.legend(['Modeled', 'Measured'])
plt.savefig(outdir + 'Modeled vs Measured PM.png', dpi = 400, bbox_inches = 'tight')

# Plot O3 model vs measured
plt.hist(o3_aqi_list, bins = 50, density = True, histtype = 'step', stacked = True, fill = False) # Modeled
plt.hist(o3_list_meas, bins = 50, density = True, histtype = 'step', stacked = True, fill = False) # Measured
plt.xlabel('AQI')
plt.ylabel('Probability')
plt.title('Measured vs Modeled Historical AQI Concentration based on O3')
plt.xlim(0, 200)
plt.legend(['Modeled', 'Measured'])
plt.savefig(outdir + 'Modeled vs Measured O3.png', dpi = 400, bbox_inches = 'tight')

# Calculate percentile values at AQI 100
dist_pm_only_meas = ECDF(pm_list_meas)
dist_pm_only_meas(100)

# Calculate AQI corresponding to the previously calculated percentile values
AQI_mod_both = np.percentile(tot_aqi_list, 99.5)
AQI_mod_PM = np.percentile(pm_aqi_list, 99.7)


# plt.legend([tot_aqi_list, pm_aqi_list, o3_aqi_list], ['Total AQI', 'PM AQI', 'O3 AQI'], loc = 'upper right')
