# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 18:23:46 2022

@author: User
"""

import geopandas as gpd
import pandas as pd
import pickle as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import numpy as np
from shapely import wkt

# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"

# Load exposure pickles
with open(indir1 + 'pm_base_REF_cell_2100', 'rb') as f:
     pm_base_REF_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_REF_cell_2050', 'rb') as f:
     pm_base_REF_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P45_cell_2100', 'rb') as f:
     pm_base_P45_cell_2100 = cp.load(f)
f.close()

# Load pickles for REF, P37, P45 adapt cases
with open(indir1 + 'pm_base_P45_cell_2050', 'rb') as f:
     pm_base_P45_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2100', 'rb') as f:
     pm_base_P37_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2050', 'rb') as f:
     pm_base_P37_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2100', 'rb') as f:
     pm_adapt_REF_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2050', 'rb') as f:
     pm_adapt_REF_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P45_cell_2100', 'rb') as f:
     pm_adapt_P45_cell_2100 = cp.load(f)
f.close()

# Load pickles for REF, P37, P45 adapt cases
with open(indir1 + 'pm_adapt_P45_cell_2050', 'rb') as f:
     pm_adapt_P45_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2100', 'rb') as f:
     pm_adapt_P37_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2050', 'rb') as f:
     pm_adapt_P37_cell_2050 = cp.load(f)
f.close()

# Load mortality pickles
with open(indir1 + 'base_burden_all', 'rb') as f:
     base_burden_all = cp.load(f)
f.close()

with open(indir1 + 'adapt_burden_all', 'rb') as f:
     adapt_burden_all = cp.load(f)
f.close()

# Load adaptation days pickles
with open(indir1 + 'adapt_count_all_P37', 'rb') as f:
     adapt_count_all_P37 = cp.load(f)
f.close()

with open(indir1 + 'adapt_count_all_P45', 'rb') as f:
     adapt_count_all_P45 = cp.load(f)
f.close()

with open(indir1 + 'adapt_count_all_REF', 'rb') as f:
     adapt_count_all_REF = cp.load(f)
f.close()

# Define parameters
# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir2 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Import population per grid cell
pop = pd.read_csv(indir2 + "benmap_pop_all.csv")
# Drop grid cells that aren't included in BenMap analysis
pop.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
pop.reset_index(inplace = True)

# Population projection factors
PPF_2050 = 1.48
PPF_2100 = 1.73

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Calculate proportion of population in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    
# Calculate weighted sums of exposure values
base_nat_2050_REF = []
base_nat_2100_REF = []   
base_nat_2050_P45 = []
base_nat_2100_P45 = []
base_nat_2050_P37 = []
base_nat_2100_P37 = []

ad_nat_2050_REF = []
ad_nat_2100_REF = []   
ad_nat_2050_P45 = []
ad_nat_2100_P45 = []
ad_nat_2050_P37 = []
ad_nat_2100_P37 = []

for r, c in rc:
    base_nat_2050_REF.append(np.mean(pm_base_REF_cell_2050[(r, c)])* pop_prop[(r, c)])
    base_nat_2100_REF.append(np.mean(pm_base_REF_cell_2100[(r, c)])* pop_prop[(r, c)])
    base_nat_2050_P45.append(np.mean(pm_base_P45_cell_2050[(r, c)])* pop_prop[(r, c)]) 
    base_nat_2100_P45.append(np.mean(pm_base_P45_cell_2100[(r, c)])* pop_prop[(r, c)]) 
    base_nat_2050_P37.append(np.mean(pm_base_P37_cell_2050[(r, c)])* pop_prop[(r, c)]) 
    base_nat_2100_P37.append(np.mean(pm_base_P37_cell_2100[(r, c)])* pop_prop[(r, c)]) 

    ad_nat_2050_REF.append(np.mean(pm_adapt_REF_cell_2050[(r, c)]) * pop_prop[(r, c)]) 
    ad_nat_2100_REF.append(np.mean(pm_adapt_REF_cell_2100[(r, c)])* pop_prop[(r, c)])     
    ad_nat_2050_P45.append(np.mean(pm_adapt_P45_cell_2050[(r, c)])* pop_prop[(r, c)])  
    ad_nat_2100_P45.append(np.mean(pm_adapt_P45_cell_2100[(r, c)])* pop_prop[(r, c)])  
    ad_nat_2050_P37.append(np.mean(pm_adapt_P37_cell_2050[(r, c)])* pop_prop[(r, c)])  
    ad_nat_2100_P37.append(np.mean(pm_adapt_P37_cell_2100[(r, c)])* pop_prop[(r, c)])  

base_nat_2050_REF = np.sum(base_nat_2050_REF)
base_nat_2100_REF = np.sum(base_nat_2100_REF)
base_nat_2050_P45 = np.sum(base_nat_2050_P45)
base_nat_2100_P45 = np.sum(base_nat_2100_P45)
base_nat_2050_P37 = np.sum(base_nat_2050_P37)
base_nat_2100_P37 = np.sum(base_nat_2100_P37)

ad_nat_2050_REF = np.sum(ad_nat_2050_REF)
ad_nat_2100_REF = np.sum(ad_nat_2100_REF)  
ad_nat_2050_P45 = np.sum(ad_nat_2050_P45)
ad_nat_2100_P45 = np.sum(ad_nat_2100_P45)
ad_nat_2050_P37 = np.sum(ad_nat_2050_P37)
ad_nat_2100_P37 = np.sum(ad_nat_2100_P37)   

# Calculate base and adaptation burden by grid cell 2050
base_burden_tot_cell_2050 = {}
adapt_burden_tot_cell_2050 = {}
base_burden_pc_cell_2050 = {}
adapt_burden_pc_cell_2050 = {}

for Pol in ['REF', 'P45', 'P37']:
    base_burden_tot_cell_2050[Pol] = {}
    adapt_burden_tot_cell_2050[Pol] = {}
    base_burden_pc_cell_2050[Pol] = {}
    adapt_burden_pc_cell_2050[Pol] = {}
    for r, c in rc:
        base_burden_tot_cell_2050[Pol][(r, c)] = []
        adapt_burden_tot_cell_2050[Pol][(r, c)] = []
        
        for Yr in year_set[0]:
            for IC in IC_list:
                base_burden_tot_cell_2050[Pol][(r, c)].append(base_burden_all[Pol][Yr][IC][(r, c)])
                adapt_burden_tot_cell_2050[Pol][(r, c)].append(adapt_burden_all[Pol][Yr][IC][(r, c)])

        base_burden_tot_cell_2050[Pol][(r, c)] = np.mean(base_burden_tot_cell_2050[Pol][(r, c)])
        adapt_burden_tot_cell_2050[Pol][(r, c)] = np.mean(adapt_burden_tot_cell_2050[Pol][(r, c)])
        base_burden_pc_cell_2050[Pol][(r, c)] = base_burden_tot_cell_2050[Pol][(r, c)] / (pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * PPF_2050)
        adapt_burden_pc_cell_2050[Pol][(r, c)] = adapt_burden_tot_cell_2050[Pol][(r, c)] / (pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * PPF_2050)

# Calculate 2050 national burden not population weighted
base_burden_tot_nat_2050 = {}
adapt_burden_tot_nat_2050 = {}
base_burden_pc_nat_2050 = {}
adapt_burden_pc_nat_2050 = {}

for Pol in ['REF', 'P45', 'P37']:
    base_burden_tot_nat_2050[Pol] = []
    adapt_burden_tot_nat_2050[Pol] = []
    base_burden_pc_nat_2050[Pol] = 0
    adapt_burden_pc_nat_2050[Pol] = 0
    for r, c in rc:
        base_burden_tot_nat_2050[Pol].append(np.mean(base_burden_tot_cell_2050[Pol][(r, c)]))
        adapt_burden_tot_nat_2050[Pol].append(np.mean(adapt_burden_tot_cell_2050[Pol][(r, c)]))
        

    base_burden_tot_nat_2050[Pol] = np.sum(base_burden_tot_nat_2050[Pol])
    adapt_burden_tot_nat_2050[Pol] = np.sum(adapt_burden_tot_nat_2050[Pol])
    base_burden_pc_nat_2050[Pol] = (base_burden_tot_nat_2050[Pol]) / (np.sum(pop['2000']) * PPF_2050)
    adapt_burden_pc_nat_2050[Pol] = (adapt_burden_tot_nat_2050[Pol]) / (np.sum(pop['2000']) * PPF_2050)
    
# Calculate 2050 national pc burden population weighted
base_burden_pc_nat_2050 = {}
adapt_burden_pc_nat_2050 = {}

for Pol in ['REF', 'P45', 'P37']:
    
    base_burden_pc_nat_2050[Pol] = []
    adapt_burden_pc_nat_2050[Pol] = []
    for r, c in rc:
        
        base_burden_pc_nat_2050[Pol].append(np.mean(base_burden_pc_cell_2050[Pol][(r, c)])* pop_prop[(r, c)])
        adapt_burden_pc_nat_2050[Pol].append(np.mean(adapt_burden_pc_cell_2050[Pol][(r, c)])* pop_prop[(r, c)])

    
    base_burden_pc_nat_2050[Pol] = np.sum(base_burden_pc_nat_2050[Pol])
    adapt_burden_pc_nat_2050[Pol] = np.sum(adapt_burden_pc_nat_2050[Pol])    
    
  
# Calculate base and adaptation burden by grid cell 2100
base_burden_tot_cell_2100 = {}
adapt_burden_tot_cell_2100 = {}
base_burden_pc_cell_2100 = {}
adapt_burden_pc_cell_2100 = {}

for Pol in ['REF', 'P45', 'P37']:
    base_burden_tot_cell_2100[Pol] = {}
    adapt_burden_tot_cell_2100[Pol] = {}
    base_burden_pc_cell_2100[Pol] = {}
    adapt_burden_pc_cell_2100[Pol] = {}
    for r, c in rc:
        base_burden_tot_cell_2100[Pol][(r, c)] = []
        adapt_burden_tot_cell_2100[Pol][(r, c)] = []
        
        for Yr in year_set[1]:
            for IC in IC_list:
                base_burden_tot_cell_2100[Pol][(r, c)].append(base_burden_all[Pol][Yr][IC][(r, c)])
                adapt_burden_tot_cell_2100[Pol][(r, c)].append(adapt_burden_all[Pol][Yr][IC][(r, c)])

        base_burden_tot_cell_2100[Pol][(r, c)] = np.mean(base_burden_tot_cell_2100[Pol][(r, c)])
        adapt_burden_tot_cell_2100[Pol][(r, c)] = np.mean(adapt_burden_tot_cell_2100[Pol][(r, c)])
        base_burden_pc_cell_2100[Pol][(r, c)] = base_burden_tot_cell_2100[Pol][(r, c)] / (pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * PPF_2100)
        adapt_burden_pc_cell_2100[Pol][(r, c)] = adapt_burden_tot_cell_2100[Pol][(r, c)] / (pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * PPF_2100)

# Calculate 2100 national burden not population weighted
base_burden_tot_nat_2100 = {}
adapt_burden_tot_nat_2100 = {}
base_burden_pc_nat_2100 = {}
adapt_burden_pc_nat_2100 = {}

for Pol in ['REF', 'P45', 'P37']:
    base_burden_tot_nat_2100[Pol] = []
    adapt_burden_tot_nat_2100[Pol] = []
    base_burden_pc_nat_2100[Pol] = 0
    adapt_burden_pc_nat_2100[Pol] = 0
    for r, c in rc:
        base_burden_tot_nat_2100[Pol].append(np.mean(base_burden_tot_cell_2100[Pol][(r, c)]))
        adapt_burden_tot_nat_2100[Pol].append(np.mean(adapt_burden_tot_cell_2100[Pol][(r, c)]))
        

    base_burden_tot_nat_2100[Pol] = np.sum(base_burden_tot_nat_2100[Pol])
    adapt_burden_tot_nat_2100[Pol] = np.sum(adapt_burden_tot_nat_2100[Pol])
    base_burden_pc_nat_2100[Pol] = (base_burden_tot_nat_2100[Pol]) / (np.sum(pop['2000']) * PPF_2100)
    adapt_burden_pc_nat_2100[Pol] = (adapt_burden_tot_nat_2100[Pol]) / (np.sum(pop['2000']) * PPF_2100)
    
# Calculate 2100 national pc burden population weighted
base_burden_pc_nat_2100 = {}
adapt_burden_pc_nat_2100 = {}

for Pol in ['REF', 'P45', 'P37']:
    
    base_burden_pc_nat_2100[Pol] = []
    adapt_burden_pc_nat_2100[Pol] = []
    for r, c in rc:
        
        base_burden_pc_nat_2100[Pol].append(np.mean(base_burden_pc_cell_2100[Pol][(r, c)])* pop_prop[(r, c)])
        adapt_burden_pc_nat_2100[Pol].append(np.mean(adapt_burden_pc_cell_2100[Pol][(r, c)])* pop_prop[(r, c)])

    
    base_burden_pc_nat_2100[Pol] = np.sum(base_burden_pc_nat_2100[Pol])
    adapt_burden_pc_nat_2100[Pol] = np.sum(adapt_burden_pc_nat_2100[Pol])    

# Calculate adaptation days per grid cell 2050
adapt_days_REF_2050_cell = {}
adapt_days_P45_2050_cell = {}
adapt_days_P37_2050_cell = {}

for r, c in rc:
    adapt_days_REF_2050_cell[(r, c)] = []
    adapt_days_P45_2050_cell[(r, c)] = []
    adapt_days_P37_2050_cell[(r, c)] = []      
    
    for Yr in year_set[0]:
        for IC in IC_list:
            adapt_days_REF_2050_cell[(r, c)].append(adapt_count_all_REF[Yr][IC][(r, c)])
            adapt_days_P45_2050_cell[(r, c)].append(adapt_count_all_P45[Yr][IC][(r, c)])
            adapt_days_P37_2050_cell[(r, c)].append(adapt_count_all_P37[Yr][IC][(r, c)])
    
    adapt_days_REF_2050_cell[(r, c)] = np.mean(adapt_days_REF_2050_cell[(r, c)])
    adapt_days_P45_2050_cell[(r, c)] = np.mean(adapt_days_P45_2050_cell)
    adapt_days_P37_2050_cell[(r, c)] = np.mean(adapt_days_P37_2050_cell)
    
adapt_days_REF_2050_tot = np.mean(adapt_days_REF_2050_cell.values())  
adapt_days_P45_2050_tot = []      
adapt_days_P37_2050_tot = []     
    
for r, c in rc:
    adapt_days_REF_2050_tot.append() 
    adapt_days_P45_2050_tot = []      
    adapt_days_P37_2050_tot = []  

