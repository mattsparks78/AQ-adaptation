# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 08:03:02 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np

# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir4 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\cb_2018_us_nation_5m\\'
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\\"

# Load Future AQ pickle
with open(indir1 + 'P37_prop_bad_2050_2100_d_pm_only', 'rb') as f:
     prop_P37 = cp.load(f)
f.close()

with open(indir1 + 'P45_prop_bad_2050_2100_d_pm_only', 'rb') as f:
     prop_P45 = cp.load(f)
f.close()

with open(indir1 + 'REF_prop_bad_2050_2100_d_pm_only', 'rb') as f:
     prop_REF = cp.load(f)
f.close()

# Load grid cells by region
with open(indir2 + 'grid_by_census_region', 'rb') as f:
     grid_region = cp.load(f)
f.close()

# Import population per grid cell
pop = pd.read_csv(indir3 + "pop_lepeule.csv")

# Create rc
rc = list(zip(pop['ROW'], pop['COL']))
    
regions = ['Northeast', 'Midwest', 'South', 'West']

reg_pop_prop = copy.deepcopy(pop)
reg_pop_prop = reg_pop_prop.drop(['2050', '2100'], axis = 1)
for reg in regions:
    reg_pop_prop[reg] = 0

# Calculate total population per region
tot_pop_reg = {}
for reg in regions:
    tot_pop_reg[reg] = []
    for r, c in rc:
        perc_in_region = grid_region.loc[(grid_region['ROW'] == r) & (grid_region['COL'] == c), reg].item()
        popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item()
        tot_pop_reg[reg].append(perc_in_region * popu)
    tot_pop_reg[reg] = sum(tot_pop_reg[reg])

# Calculate population proportion for cells per region
for reg in regions:
    tot_pop = tot_pop_reg[reg]
    reg_pop_prop[reg] = reg_pop_prop['2000'] * grid_region[reg] / tot_pop
        
# Caclulate time series of regional average bad air days REF 2050
# Create one dict to hold them all
reg_prop_REF_2050 = {}  

for reg in regions:
    reg_prop_REF_2050[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_REF_2050[reg][IC] = []
        for Yr in range(2036,2066):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_REF[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_REF_2050[reg][IC].append(sum(a))
title = 'High AQI days Per Region (REF 2050 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2050['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_REF_2050['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_REF_2050['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_REF_2050['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
# fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Caclulate time series of regional average bad air days P45 2050
# Create one dict to hold them all
reg_prop_P45_2050 = {}  

for reg in regions:
    reg_prop_P45_2050[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_P45_2050[reg][IC] = []
        for Yr in range(2036,2066):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_P45[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_P45_2050[reg][IC].append(sum(a))
title = 'High AQI days Per Region (P45 2050 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_P45_2050['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_P45_2050['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_P45_2050['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_P45_2050['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
# fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Caclulate time series of regional average bad air days P37 2050
# Create one dict to hold them all
reg_prop_P37_2050 = {}  

for reg in regions:
    reg_prop_P37_2050[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_P37_2050[reg][IC] = []
        for Yr in range(2036,2066):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_P37[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_P37_2050[reg][IC].append(sum(a))
title = 'High AQI days Per Region (P37 2050 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_P37_2050['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_P37_2050['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_P37_2050['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_P37_2050['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
# fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Check mean and SD across IC, years

pREF = {}
pP45 = {}
pP37 = {}

for reg in regions:
    pREF[reg] = []
    pP45[reg] = []
    pP37[reg] = []
    for IC in ['IC' + str(i) for i in range(1,6)]:
        pREF[reg].extend(reg_prop_REF_2050[reg][IC])
        pP45[reg].extend(reg_prop_P45_2050[reg][IC])
        pP37[reg].extend(reg_prop_P37_2050[reg][IC])
    
    pREF[reg] = np.array(pREF[reg])
    pP45[reg] = np.array(pP45[reg])
    pP37[reg] = np.array(pP37[reg])

meanREF = {}
meanP45 = {}
meanP37 = {}
sdREF = {}
sdP45 = {}
sdP37 = {}

for reg in regions:
    meanREF[reg] = np.mean(pREF[reg])
    meanP45[reg] = np.mean(pP45[reg])
    meanP37[reg] = np.mean(pP37[reg])
    
    sdREF[reg] = np.std(pREF[reg])
    sdP45[reg] = np.std(pP45[reg])
    sdP37[reg] = np.std(pP37[reg])


# Caclulate time series of regional average bad air days REF 2100
# Create one dict to hold them all
reg_prop_REF_2100 = {}  

for reg in regions:
    reg_prop_REF_2100[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_REF_2100[reg][IC] = []
        for Yr in range(2086,2116):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_REF[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_REF_2100[reg][IC].append(sum(a))
title = 'High AQI days Per Region (REF 2100 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2100['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_REF_2100['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_REF_2100['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_REF_2100['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Caclulate time series of regional average bad air days P45 2100
# Create one dict to hold them all
reg_prop_P45_2100 = {}  

for reg in regions:
    reg_prop_P45_2100[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_P45_2100[reg][IC] = []
        for Yr in range(2086,2116):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_P45[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_P45_2100[reg][IC].append(sum(a))
title = 'High AQI days Per Region (REF 2100 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_P45_2100['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_P45_2100['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_P45_2100['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_P45_2100['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Caclulate time series of regional average bad air days P37 2100
# Create one dict to hold them all
reg_prop_P37_2100 = {}  

for reg in regions:
    reg_prop_P37_2100[reg] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        reg_prop_P37_2100[reg][IC] = []
        for Yr in range(2086,2116):
            a = []
            for r, c in rc:
                popu_prop = reg_pop_prop.loc[(reg_pop_prop['ROW'] == r) & (reg_pop_prop['COL'] == c), reg].item()
                aqi_prop = prop_P37[str(Yr)][IC][(r, c)]
                a.append(popu_prop * aqi_prop)
            reg_prop_P37_2100[reg][IC].append(sum(a))
title = 'High AQI days Per Region (REF 2100 IC1)'

fig, ax = plt.subplots()
plt.plot(reg_prop_P37_2100['Northeast']['IC1'], color = 'green', label = 'Northeast')
plt.plot(reg_prop_P37_2100['Midwest']['IC1'], color = 'red', label = 'Midwest')
plt.plot(reg_prop_P37_2100['South']['IC1'], color = 'blue', label = 'South')
plt.plot(reg_prop_P37_2100['West']['IC1'], color = 'black', label = 'West')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Plot all scenarios for Northeast region
title = 'High AQI Days in Northeast Region by Scenario'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2050['Northeast']['IC1'], color = 'green', label = 'REF 2050')
plt.plot(reg_prop_P45_2050['Northeast']['IC1'], color = 'red', label = 'P45 2050')
plt.plot(reg_prop_P37_2050['Northeast']['IC1'], color = 'blue', label = 'P37 2050')
plt.plot(reg_prop_REF_2100['Northeast']['IC1'], color = 'black', label = 'REF 2100')
plt.plot(reg_prop_P45_2100['Northeast']['IC1'], color = 'yellow', label = 'P45 2100')
plt.plot(reg_prop_P37_2100['Northeast']['IC1'], color = 'purple', label = 'P37 2100')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Plot all scenarios for Midwest region
title = 'High AQI Days in Midwest Region by Scenario'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2050['Midwest']['IC1'], color = 'green', label = 'REF 2050')
plt.plot(reg_prop_P45_2050['Midwest']['IC1'], color = 'red', label = 'P45 2050')
plt.plot(reg_prop_P37_2050['Midwest']['IC1'], color = 'blue', label = 'P37 2050')
plt.plot(reg_prop_REF_2100['Midwest']['IC1'], color = 'black', label = 'REF 2100')
plt.plot(reg_prop_P45_2100['Midwest']['IC1'], color = 'yellow', label = 'P45 2100')
plt.plot(reg_prop_P37_2100['Midwest']['IC1'], color = 'purple', label = 'P37 2100')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Plot all scenarios for South region
title = 'High AQI Days in South Region by Scenario'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2050['South']['IC1'], color = 'green', label = 'REF 2050')
plt.plot(reg_prop_P45_2050['South']['IC1'], color = 'red', label = 'P45 2050')
plt.plot(reg_prop_P37_2050['South']['IC1'], color = 'blue', label = 'P37 2050')
plt.plot(reg_prop_REF_2100['South']['IC1'], color = 'black', label = 'REF 2100')
plt.plot(reg_prop_P45_2100['South']['IC1'], color = 'yellow', label = 'P45 2100')
plt.plot(reg_prop_P37_2100['South']['IC1'], color = 'purple', label = 'P37 2100')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

# Plot all scenarios for West region
title = 'High AQI Days in West Region by Scenario'

fig, ax = plt.subplots()
plt.plot(reg_prop_REF_2050['West']['IC1'], color = 'green', label = 'REF 2050')
plt.plot(reg_prop_P45_2050['West']['IC1'], color = 'red', label = 'P45 2050')
plt.plot(reg_prop_P37_2050['West']['IC1'], color = 'blue', label = 'P37 2050')
plt.plot(reg_prop_REF_2100['West']['IC1'], color = 'black', label = 'REF 2100')
plt.plot(reg_prop_P45_2100['West']['IC1'], color = 'yellow', label = 'P45 2100')
plt.plot(reg_prop_P37_2100['West']['IC1'], color = 'purple', label = 'P37 2100')
ax.set_xlabel('Simulation Year')
ax.set_ylabel('Pop-Weighted Proportion of Days AQI > 100 (PM Only')
ax.set_title(title)
ax.legend()
fig.savefig(outdir + title + '.png', dpi =400, bbox_inches = 'tight')

