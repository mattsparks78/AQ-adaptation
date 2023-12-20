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
indir1 = r".\Data\Rational Actor Infiltration Dem\\"
indir2 = r".\Data\General\\"
indir3 = r".\Data\AGU2020-Fig3\\"
outdir1 = r".\Data\Valuations All Infiltration Dem\\"
outdir2 = r".\Data\Rational Actor Infiltration Dem\\"

# Set Infiltration Rate
inf_rate = 0.2

# Load exposure pickles
with open(indir1 + 'pm_base_REF_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_base_REF_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_REF_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_base_REF_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P45_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_base_P45_cell_2100 = cp.load(f)
f.close()

# Load pickles for REF, P37, P45 adapt cases
with open(indir1 + 'pm_base_P45_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_base_P45_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_base_P37_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_base_P37_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_base_P37_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_adapt_REF_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_REF_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_adapt_REF_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P45_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_adapt_P45_cell_2100 = cp.load(f)
f.close()

# Load pickles for REF, P37, P45 adapt cases
with open(indir1 + 'pm_adapt_P45_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_adapt_P45_cell_2050 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2100_' + str(inf_rate), 'rb') as f:
     pm_adapt_P37_cell_2100 = cp.load(f)
f.close()

with open(indir1 + 'pm_adapt_P37_cell_2050_' + str(inf_rate), 'rb') as f:
     pm_adapt_P37_cell_2050 = cp.load(f)
f.close()

# Load mortality pickles
with open(indir1 + 'base_burden_all_' + str(inf_rate) , 'rb') as f:
     base_burden_all = cp.load(f)
f.close()

with open(indir1 + 'adapt_burden_all_' + str(inf_rate), 'rb') as f:
     adapt_burden_all = cp.load(f)
f.close()

# Load adaptation days pickles
with open(indir1 + 'adapt_count_all_P37_' + str(inf_rate), 'rb') as f:
     adapt_count_all_P37 = cp.load(f)
f.close()

with open(indir1 + 'adapt_count_all_P45_' + str(inf_rate), 'rb') as f:
     adapt_count_all_P45 = cp.load(f)
f.close()

with open(indir1 + 'adapt_count_all_REF_' + str(inf_rate), 'rb') as f:
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

# Population
pop = pd.read_csv(indir2 + 'pop_lepeule_race.csv')

# List of cells
rc = list(zip(pop['ROW'], pop['COL']))

# Population projection factors
PPF_2050 = 1.48
PPF_2100 = 1.73

# Calculate proportion of population in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    
# Create race dictionary
grid_race_dict = {}
for r, c in rc:
    grid_race_dict[(r, c)] = {}
    tot_list = []
    for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        grid_race_dict[(r, c)][race] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), race].item()
        tot_list.append(pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), race].item())
    grid_race_dict[(r, c)]['Total'] = sum(tot_list)
    
tot_pop_list = []
for r, c in rc:
    tot_pop_list.append(grid_race_dict[(r, c)]['Total'])
tot_pop = sum(tot_pop_list)

# Create dict with grid/race breakdown
tot_pop_list = []
for r, c in rc:
    tot_pop_list.append(grid_race_dict[(r, c)]['Total'])
tot_pop = sum(tot_pop_list)

# Combine others into all others
for r, c in rc:
    grid_race_dict[(r, c)]['All Others'] = grid_race_dict[(r, c)]['Asian'] + grid_race_dict[(r, c)]['Native American'] + \
        grid_race_dict[(r, c)]['Pacific'] + grid_race_dict[(r, c)]['Other'] + grid_race_dict[(r, c)]['Two Plus']
        
tot_race_pop = {}
for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    tot_race_pop[race] = []
    for r, c in rc:
        tot_race_pop[race].append(grid_race_dict[(r, c)][race])
    tot_race_pop[race] = sum(tot_race_pop[race])
    
grid_race_prop = {}
for r, c in rc:
    grid_race_prop[(r, c)] = {}
    for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        grid_race_prop[(r, c)][race] = grid_race_dict[(r, c)][race] / tot_race_pop[race]   

#%%
# Calculate weighted sums of exposure values
base_nat_2050_REF = []
base_nat_2100_REF = []   
base_nat_2050_P45 = []
base_nat_2100_P45 = []
base_nat_2050_P37 = []
base_nat_2100_P37 = []

for r, c in rc:
    base_nat_2050_REF.append(np.mean(pm_base_REF_cell_2050['White'][(r, c)]) * grid_race_prop[(r, c)][race])
    base_nat_2100_REF.append(np.mean(pm_base_REF_cell_2100['White'][(r, c)]) * grid_race_prop[(r, c)][race])
    base_nat_2050_P45.append(np.mean(pm_base_P45_cell_2050['White'][(r, c)]) * grid_race_prop[(r, c)][race]) 
    base_nat_2100_P45.append(np.mean(pm_base_P45_cell_2100['White'][(r, c)]) * grid_race_prop[(r, c)][race]) 
    base_nat_2050_P37.append(np.mean(pm_base_P37_cell_2050['White'][(r, c)]) * grid_race_prop[(r, c)][race]) 
    base_nat_2100_P37.append(np.mean(pm_base_P37_cell_2100['White'][(r, c)]) * grid_race_prop[(r, c)][race]) 
    
base_nat_2050_REF = np.sum(base_nat_2050_REF)
base_nat_2100_REF = np.sum(base_nat_2100_REF)
base_nat_2050_P45 = np.sum(base_nat_2050_P45)
base_nat_2100_P45 = np.sum(base_nat_2100_P45)
base_nat_2050_P37 = np.sum(base_nat_2050_P37)
base_nat_2100_P37 = np.sum(base_nat_2100_P37)

    
ad_nat_2050_REF = {}
ad_nat_2100_REF = {}   
ad_nat_2050_P45 = {}
ad_nat_2100_P45 = {}
ad_nat_2050_P37 = {}
ad_nat_2100_P37 = {}

for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    ad_nat_2050_REF[race] = []
    ad_nat_2100_REF[race] = []  
    ad_nat_2050_P45[race] = []
    ad_nat_2100_P45[race] = []
    ad_nat_2050_P37[race] = []
    ad_nat_2100_P37[race] = []
    for r, c in rc:
        ad_nat_2050_REF[race].append(np.mean(pm_adapt_REF_cell_2050[race][(r, c)]) * grid_race_prop[(r, c)][race]) 
        ad_nat_2100_REF[race].append(np.mean(pm_adapt_REF_cell_2100[race][(r, c)])* grid_race_prop[(r, c)][race])     
        ad_nat_2050_P45[race].append(np.mean(pm_adapt_P45_cell_2050[race][(r, c)])* grid_race_prop[(r, c)][race])  
        ad_nat_2100_P45[race].append(np.mean(pm_adapt_P45_cell_2100[race][(r, c)])* grid_race_prop[(r, c)][race])  
        ad_nat_2050_P37[race].append(np.mean(pm_adapt_P37_cell_2050[race][(r, c)])* grid_race_prop[(r, c)][race])  
        ad_nat_2100_P37[race].append(np.mean(pm_adapt_P37_cell_2100[race][(r, c)])* grid_race_prop[(r, c)][race])
            
    ad_nat_2050_REF[race] = np.sum(ad_nat_2050_REF[race])
    ad_nat_2100_REF[race] = np.sum(ad_nat_2100_REF[race])  
    ad_nat_2050_P45[race] = np.sum(ad_nat_2050_P45[race])
    ad_nat_2100_P45[race] = np.sum(ad_nat_2100_P45[race])
    ad_nat_2050_P37[race] = np.sum(ad_nat_2050_P37[race])
    ad_nat_2100_P37[race] = np.sum(ad_nat_2100_P37[race])                 


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
        base_burden_tot_cell_2050[Pol][(r, c)] = {}
        adapt_burden_tot_cell_2050[Pol][(r, c)] = {}
        base_burden_pc_cell_2050[Pol][(r, c)] = {}
        adapt_burden_pc_cell_2050[Pol][(r, c)] = {}
        for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
            base_burden_tot_cell_2050[Pol][(r, c)][race] = []
            adapt_burden_tot_cell_2050[Pol][(r, c)][race] = []
            
            for Yr in year_set[0]:
                for IC in IC_list:
                    base_burden_tot_cell_2050[Pol][(r, c)][race].append(base_burden_all[Pol][Yr][IC][(r, c)][race])
                    adapt_burden_tot_cell_2050[Pol][(r, c)][race].append(adapt_burden_all[Pol][Yr][IC][(r, c)][race])
    
            base_burden_tot_cell_2050[Pol][(r, c)][race] = np.mean(base_burden_tot_cell_2050[Pol][(r, c)][race])
            adapt_burden_tot_cell_2050[Pol][(r, c)][race] = np.mean(adapt_burden_tot_cell_2050[Pol][(r, c)][race])
            base_burden_pc_cell_2050[Pol][(r, c)][race] = base_burden_tot_cell_2050[Pol][(r, c)][race] / (grid_race_dict[(r, c)][race] * PPF_2050)
            adapt_burden_pc_cell_2050[Pol][(r, c)][race] = adapt_burden_tot_cell_2050[Pol][(r, c)][race] / (grid_race_dict[(r, c)][race] * PPF_2050)

    
# Calculate 2050 national pc burden population weighted
base_burden_pc_nat_2050 = {}
adapt_burden_pc_nat_2050 = {}

for Pol in ['REF', 'P45', 'P37']:
    base_burden_pc_nat_2050[Pol] = {}
    adapt_burden_pc_nat_2050[Pol] = {}
    for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        base_burden_pc_nat_2050[Pol][race] = []
        adapt_burden_pc_nat_2050[Pol][race] = []
        for r, c in rc:
            
            base_burden_pc_nat_2050[Pol][race].append(np.mean(base_burden_pc_cell_2050[Pol][(r, c)][race])* grid_race_prop[(r, c)][race])
            adapt_burden_pc_nat_2050[Pol][race].append(np.mean(adapt_burden_pc_cell_2050[Pol][(r, c)][race])* grid_race_prop[(r, c)][race])
    
        
        base_burden_pc_nat_2050[Pol][race] = np.sum(base_burden_pc_nat_2050[Pol][race])
        adapt_burden_pc_nat_2050[Pol][race] = np.sum(adapt_burden_pc_nat_2050[Pol][race])    
    
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
        base_burden_tot_cell_2100[Pol][(r, c)] = {}
        adapt_burden_tot_cell_2100[Pol][(r, c)] = {}
        base_burden_pc_cell_2100[Pol][(r, c)] = {}
        adapt_burden_pc_cell_2100[Pol][(r, c)] = {}
        for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
            base_burden_tot_cell_2100[Pol][(r, c)][race] = []
            adapt_burden_tot_cell_2100[Pol][(r, c)][race] = []
            for Yr in year_set[1]:
                for IC in IC_list:
                    base_burden_tot_cell_2100[Pol][(r, c)][race].append(base_burden_all[Pol][Yr][IC][(r, c)][race])
                    adapt_burden_tot_cell_2100[Pol][(r, c)][race].append(adapt_burden_all[Pol][Yr][IC][(r, c)][race])
    
            base_burden_tot_cell_2100[Pol][(r, c)][race] = np.mean(base_burden_tot_cell_2100[Pol][(r, c)][race])
            adapt_burden_tot_cell_2100[Pol][(r, c)][race] = np.mean(adapt_burden_tot_cell_2100[Pol][(r, c)][race])
            base_burden_pc_cell_2100[Pol][(r, c)][race] = base_burden_tot_cell_2100[Pol][(r, c)][race] / (grid_race_dict[(r, c)][race] * PPF_2100)
            adapt_burden_pc_cell_2100[Pol][(r, c)][race] = adapt_burden_tot_cell_2100[Pol][(r, c)][race] / (grid_race_dict[(r, c)][race] * PPF_2100)

# Calculate 2100 national pc burden population weighted
base_burden_pc_nat_2100 = {}
adapt_burden_pc_nat_2100 = {}
for Pol in ['REF', 'P45', 'P37']:
    base_burden_pc_nat_2100[Pol] = {}
    adapt_burden_pc_nat_2100[Pol] = {}
    for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        base_burden_pc_nat_2100[Pol][race] = []
        adapt_burden_pc_nat_2100[Pol][race] = []
        for r, c in rc:
            base_burden_pc_nat_2100[Pol][race].append(np.mean(base_burden_pc_cell_2100[Pol][(r, c)][race])* grid_race_prop[(r, c)][race])
            adapt_burden_pc_nat_2100[Pol][race].append(np.mean(adapt_burden_pc_cell_2100[Pol][(r, c)][race])* grid_race_prop[(r, c)][race])
        base_burden_pc_nat_2100[Pol][race] = np.sum(base_burden_pc_nat_2100[Pol][race])
        adapt_burden_pc_nat_2100[Pol][race] = np.sum(adapt_burden_pc_nat_2100[Pol][race])    

# #Output pickles
# with open(outdir1 + 'base_burden_pc_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_pc_cell_2050,f)
# f.close()

# with open(outdir1 + 'base_burden_pc_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_pc_cell_2100,f)
# f.close()

# with open(outdir1 + 'base_burden_pc_nat_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_pc_nat_2050,f)
# f.close()

# with open(outdir1 + 'base_burden_pc_nat_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_pc_nat_2100,f)
# f.close()

# with open(outdir1 + 'base_burden_tot_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_tot_cell_2050,f)
# f.close()

# with open(outdir1 + 'base_burden_tot_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden_tot_cell_2100,f)
# f.close()

# # Output adaptation pickles
# with open(outdir2 + 'adapt_burden_pc_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_pc_cell_2050,f)
# f.close()

# with open(outdir2 + 'adapt_burden_pc_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_pc_cell_2100,f)
# f.close()

# with open(outdir2 + 'adapt_burden_pc_nat_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_pc_nat_2050,f)
# f.close()

# with open(outdir2 + 'adapt_burden_pc_nat_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_pc_nat_2100,f)
# f.close()

# with open(outdir2 + 'adapt_burden_tot_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_tot_cell_2050,f)
# f.close()

# with open(outdir2 + 'adapt_burden_tot_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden_tot_cell_2100,f)
# f.close()







