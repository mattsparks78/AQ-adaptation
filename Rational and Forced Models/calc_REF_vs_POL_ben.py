# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 06:22:34 2022

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
indir1 = r".\Data\General\\"
indir2 = r".\Data\Rational Actor Infiltration\\"
indir3 = r".\Data\AGU2020-Fig3\\"
outdir = r".\Data\Rational Actor Infiltration\\"

# Set Infiltration Rate
inf_rate = 0.2


# Load pickles for REF, P37, P45 baseline cases
with open(indir2 + 'PM_out_base_REF_' + str(inf_rate), 'rb') as f:
     pm_base_REF = cp.load(f)
f.close()

with open(indir2 + 'PM_out_base_P45_' + str(inf_rate), 'rb') as f:
     pm_base_P45 = cp.load(f)
f.close()

with open(indir2 + 'PM_out_base_P37_' + str(inf_rate), 'rb') as f:
     pm_base_P37 = cp.load(f)
f.close()

# Load pickles for REF, P37, P45 adapt cases
with open(indir2 + 'PM_out_adapt_REF_' + str(inf_rate), 'rb') as f:
     pm_adapt_REF = cp.load(f)
f.close()

with open(indir2 + 'PM_out_adapt_P45_' + str(inf_rate), 'rb') as f:
     pm_adapt_P45 = cp.load(f)
f.close()

with open(indir2 + 'PM_out_adapt_P37_' + str(inf_rate), 'rb') as f:
     pm_adapt_P37 = cp.load(f)
f.close()

# Define parameters
# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis = 0, inplace = True)
# Reset index back to starting at 0
wage.reset_index(inplace = True)

# Import population per grid cell
pop = pd.read_csv(indir1 + "pop_lepeule.csv")

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']),list(wage['COL'])))

# Set adaptation time in hours. Assuming full adaptation because it maximizes utility.
AT = 1

#Set relative risk
RR = 1.14
# Calculate attributable fraction 
AF = (RR - 1)/ RR

# Set VSL. Converted from $2006 USD to 2005
VSL = 7104000

# Create dict to hold health burden data
base_burden = {}
adapt_burden = {}
dif_burden = {}

# Calculate burdens for each case
for Pol in ['REF', 'P37', 'P45']:
    base_burden[Pol] = {}
    adapt_burden[Pol] = {}
    dif_burden[Pol] = {}
    for years in year_set:
        # Iterate over all years in year set
        if np.mean(years) < 2070:
            if Pol == 'REF':
                # Economic projection factor
                EPF = 2.01
            elif Pol == 'P37':
                EPF = 1.82
            elif Pol == 'P45':
                EPF = 1.86
            # Population projection factor
            PPF = 1.48
            # PM Mortality Factor (for baseline mortaltiy)
            PMMF = 1.12
            # Set VSL. Converted from $2006 USD to 2005
            VSL = 7104000
            VSL = VSL * EPF
        else:
            if Pol == 'REF':
                EPF = 3.94
            elif Pol == 'P37':
                EPF = 3.37
            elif Pol == 'P45':
                EPF = 3.5
            PPF = 1.73
            PMMF = 0.92
            # Set VSL. Converted from $2006 USD to 2005
            VSL = 7104000
            VSL = VSL * EPF
        for Yr in years:
            # Initialize next level of dictionaries
            base_burden[Pol][Yr] = {}
            adapt_burden[Pol][Yr] = {}
            dif_burden[Pol][Yr] = {}
            for IC in IC_list:
                base_burden[Pol][Yr][IC] = {}
                adapt_burden[Pol][Yr][IC] = {}
                dif_burden[Pol][Yr][IC] = {}
                for r, c in rc:
                    # Calculate burdens
                    bmr = bm.loc[(bm['ROW'] == r) & (bm['COL'] ==c), 'BaseMortRate'].item() * PMMF
                    popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() * PPF
                    
                    # Select correct scenario
                    if Pol == 'REF':
                        dPM_base = pm_base_REF[Yr][IC][(r, c)]
                        dPM_adapt = pm_adapt_REF[Yr][IC][(r, c)]
                    elif Pol == 'P37':
                        dPM_base = pm_base_P37[Yr][IC][(r, c)]
                        dPM_adapt = pm_adapt_P37[Yr][IC][(r, c)]
                    elif Pol == 'P45':
                        dPM_base = pm_base_P45[Yr][IC][(r, c)]
                        dPM_adapt = pm_adapt_P45[Yr][IC][(r, c)]

                    base_burden[Pol][Yr][IC][(r, c)] = bmr * AF * dPM_base * 0.1 * popu * VSL
                    adapt_burden[Pol][Yr][IC][(r, c)] = bmr * AF * dPM_adapt * 0.1 * popu * VSL

                    dif_burden[Pol][Yr][IC][(r, c)] = base_burden[Pol][Yr][IC][(r, c)] - adapt_burden[Pol][Yr][IC][(r, c)]


# Calculate P45 policy cobenefits 
coben_P45 = {}
adapt_coben_P45 = {}

for years in year_set:
    for Yr in years:
        coben_P45[Yr] = {}
        adapt_coben_P45[Yr] = {}
        for IC in IC_list:
            coben_P45[Yr][IC] = {}
            adapt_coben_P45[Yr][IC] = {}
            for r, c in rc:
                coben_P45[Yr][IC][(r, c)] = base_burden['REF'][Yr][IC][(r, c)] - base_burden['P45'][Yr][IC][(r, c)]
                adapt_coben_P45[Yr][IC][(r, c)] = adapt_burden['REF'][Yr][IC][(r, c)] - adapt_burden['P45'][Yr][IC][(r, c)]
               
# Calculate P37 policy cobenefits 
coben_P37 = {}
adapt_coben_P37 = {}

for years in year_set:
    for Yr in years:
        coben_P37[Yr] = {}
        adapt_coben_P37[Yr] = {}
        for IC in IC_list:
            coben_P37[Yr][IC] = {}
            adapt_coben_P37[Yr][IC] = {}
            for r, c in rc:
                coben_P37[Yr][IC][(r, c)] = base_burden['REF'][Yr][IC][(r, c)] - base_burden['P37'][Yr][IC][(r, c)]
                adapt_coben_P37[Yr][IC][(r, c)] = adapt_burden['REF'][Yr][IC][(r, c)] - adapt_burden['P37'][Yr][IC][(r, c)]
                
# Calculate P45 policy cobenefits by cell for 2050
coben_P45_2050_cell = {}
adapt_coben_P45_2050_cell = {}
dif_coben_cell_P45_2050 = {}
dif_coben_mean_cell_P45_2050 = {}
for r, c in rc:               
    coben = []
    ad_coben = []
    for Yr in year_set[0]:
        for IC in IC_list:
            coben.append(coben_P45[Yr][IC][(r, c)])
            ad_coben.append(adapt_coben_P45[Yr][IC][(r, c)])
    coben_P45_2050_cell[(r, c)] = np.array(coben)
    adapt_coben_P45_2050_cell[(r, c)] = np.array(ad_coben)                    
    dif_coben_cell_P45_2050[(r, c)] = (coben_P45_2050_cell[(r, c)] - adapt_coben_P45_2050_cell[(r, c)]) / coben_P45_2050_cell[(r, c)] 
    dif_coben_mean_cell_P45_2050[(r, c)] = np.mean(dif_coben_cell_P45_2050[(r, c)])

# Calculate P45 policy cobenefits by cell for 2100
coben_P45_2100_cell = {}
adapt_coben_P45_2100_cell = {}
dif_coben_cell_P45_2100 = {}
dif_coben_mean_cell_P45_2100 = {}
for r, c in rc:               
    coben = []
    ad_coben = []
    for Yr in year_set[1]:
        for IC in IC_list:
            coben.append(coben_P45[Yr][IC][(r, c)])
            ad_coben.append(adapt_coben_P45[Yr][IC][(r, c)])
    coben_P45_2100_cell[(r, c)] = np.array(coben)
    adapt_coben_P45_2100_cell[(r, c)] = np.array(ad_coben)                    
    dif_coben_cell_P45_2100[(r, c)] = (coben_P45_2100_cell[(r, c)] - adapt_coben_P45_2100_cell[(r, c)]) / coben_P45_2100_cell[(r, c)] 
    dif_coben_mean_cell_P45_2100[(r, c)] = np.mean(dif_coben_cell_P45_2100[(r, c)])
    
# Calculate P37 policy cobenefits by cell for 2050
coben_P37_2050_cell = {}
adapt_coben_P37_2050_cell = {}
dif_coben_cell_P37_2050 = {}
dif_coben_mean_cell_P37_2050 = {}
for r, c in rc:               
    coben = []
    ad_coben = []
    for Yr in year_set[0]:
        for IC in IC_list:
            coben.append(coben_P37[Yr][IC][(r, c)])
            ad_coben.append(adapt_coben_P37[Yr][IC][(r, c)])
    coben_P37_2050_cell[(r, c)] = np.array(coben)
    adapt_coben_P37_2050_cell[(r, c)] = np.array(ad_coben)                    
    dif_coben_cell_P37_2050[(r, c)] = (coben_P37_2050_cell[(r, c)] - adapt_coben_P37_2050_cell[(r, c)]) / coben_P37_2050_cell[(r, c)] 
    dif_coben_mean_cell_P37_2050[(r, c)] = np.mean(dif_coben_cell_P37_2050[(r, c)])
    
# Calculate P37 policy cobenefits by cell for 2100
coben_P37_2100_cell = {}
adapt_coben_P37_2100_cell = {}
dif_coben_cell_P37_2100 = {}
dif_coben_mean_cell_P37_2100 = {}
for r, c in rc:               
    coben = []
    ad_coben = []
    for Yr in year_set[1]:
        for IC in IC_list:
            coben.append(coben_P37[Yr][IC][(r, c)])
            ad_coben.append(adapt_coben_P37[Yr][IC][(r, c)])
    coben_P37_2100_cell[(r, c)] = np.array(coben)
    adapt_coben_P37_2100_cell[(r, c)] = np.array(ad_coben)                    
    dif_coben_cell_P37_2100[(r, c)] = (coben_P37_2100_cell[(r, c)] - adapt_coben_P37_2100_cell[(r, c)]) / coben_P37_2100_cell[(r, c)] 
    dif_coben_mean_cell_P37_2100[(r, c)] = np.mean(dif_coben_cell_P37_2100[(r, c)])

# Calculate national cobenefits
nat_coben = {}

for year in [2050, 2100]:
    nat_coben[year] = {}
    for pol in ['P37', 'P45']:
        if year == 2050:
            if pol == 'P37':
                cell_dict = coben_P37_2050_cell
            elif pol == 'P45':
                cell_dict = coben_P45_2050_cell
        if year == 2100:
            if pol == 'P37':
                cell_dict = coben_P37_2100_cell
            elif pol == 'P45':
                cell_dict = coben_P45_2100_cell
                
        for r, c in rc:
            cell_dict[(r, c)] = np.mean(cell_dict[(r, c)])
        nat_coben[year][pol] = sum(cell_dict.values())* 1.4994 / 1E12




    
# Calculate base and adaptation pm by grid cell
# Create dictionaries
pm_adapt_P37_cell_2050 = {}
pm_adapt_P45_cell_2050 = {}
pm_adapt_REF_cell_2050 = {}
pm_base_P37_cell_2050 = {}
pm_base_P45_cell_2050 = {}
pm_base_REF_cell_2050 = {}

# Iterate over grid cells for 2050
for r, c in rc:
    ad_P37_cell_2050 = []
    ad_P45_cell_2050 = []
    ad_REF_cell_2050 = []
    b_P37_cell_2050 = []
    b_P45_cell_2050 = []
    b_REF_cell_2050 = []
    for Yr in year_set[0]:
        for IC in IC_list:
            ad_P37_cell_2050.append(pm_adapt_P37[Yr][IC][(r, c)])
            ad_P45_cell_2050.append(pm_adapt_P45[Yr][IC][(r, c)])
            ad_REF_cell_2050.append(pm_adapt_REF[Yr][IC][(r, c)])
            b_P37_cell_2050.append(pm_base_P37[Yr][IC][(r, c)])
            b_P45_cell_2050.append(pm_base_P45[Yr][IC][(r, c)])
            b_REF_cell_2050.append(pm_base_REF[Yr][IC][(r, c)])

    pm_adapt_P37_cell_2050[(r, c)] = np.array(ad_P37_cell_2050)
    pm_adapt_P45_cell_2050[(r, c)] = np.array(ad_P45_cell_2050)
    pm_adapt_REF_cell_2050[(r, c)] = np.array(ad_REF_cell_2050)
    pm_base_P37_cell_2050[(r, c)] = np.array(b_P37_cell_2050)
    pm_base_P45_cell_2050[(r, c)] = np.array(b_P45_cell_2050)
    pm_base_REF_cell_2050[(r, c)] = np.array(b_REF_cell_2050)
    
# Do the same for 2100
# Create dictionaries
pm_adapt_P37_cell_2100 = {}
pm_adapt_P45_cell_2100 = {}
pm_adapt_REF_cell_2100 = {}
pm_base_P37_cell_2100 = {}
pm_base_P45_cell_2100 = {}
pm_base_REF_cell_2100 = {}

# Iterate over grid cells for 2100
for r, c in rc:
    ad_P37_cell_2100 = []
    ad_P45_cell_2100 = []
    ad_REF_cell_2100 = []
    b_P37_cell_2100 = []
    b_P45_cell_2100 = []
    b_REF_cell_2100 = []
    for Yr in year_set[1]:
        for IC in IC_list:
            ad_P37_cell_2100.append(pm_adapt_P37[Yr][IC][(r, c)])
            ad_P45_cell_2100.append(pm_adapt_P45[Yr][IC][(r, c)])
            ad_REF_cell_2100.append(pm_adapt_REF[Yr][IC][(r, c)])
            b_P37_cell_2100.append(pm_base_P37[Yr][IC][(r, c)])
            b_P45_cell_2100.append(pm_base_P45[Yr][IC][(r, c)])
            b_REF_cell_2100.append(pm_base_REF[Yr][IC][(r, c)])

    pm_adapt_P37_cell_2100[(r, c)] = np.array(ad_P37_cell_2100)
    pm_adapt_P45_cell_2100[(r, c)] = np.array(ad_P45_cell_2100)
    pm_adapt_REF_cell_2100[(r, c)] = np.array(ad_REF_cell_2100)
    pm_base_P37_cell_2100[(r, c)] = np.array(b_P37_cell_2100)
    pm_base_P45_cell_2100[(r, c)] = np.array(b_P45_cell_2100)
    pm_base_REF_cell_2100[(r, c)] = np.array(b_REF_cell_2100)
    
    
# Output concentration pickles
# 2050
# with open(outdir + 'pm_adapt_P37_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_P37_cell_2050,f)
# f.close()

# with open(outdir + 'pm_adapt_P45_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_P45_cell_2050,f)
# f.close()

# with open(outdir + 'pm_adapt_REF_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_REF_cell_2050,f)
# f.close()

# with open(outdir + 'pm_base_P37_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_P37_cell_2050,f)
# f.close()

# with open(outdir + 'pm_base_P45_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_P45_cell_2050,f)
# f.close()

# with open(outdir + 'pm_base_REF_cell_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_REF_cell_2050,f)
# f.close()

# # # 2100
# with open(outdir + 'pm_adapt_P37_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_P37_cell_2100,f)
# f.close()

# with open(outdir + 'pm_adapt_P45_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_P45_cell_2100,f)
# f.close()

# with open(outdir + 'pm_adapt_REF_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_adapt_REF_cell_2100,f)
# f.close()

# with open(outdir + 'pm_base_P37_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_P37_cell_2100,f)
# f.close()

# with open(outdir + 'pm_base_P45_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_P45_cell_2100,f)
# f.close()

# with open(outdir + 'pm_base_REF_cell_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(pm_base_REF_cell_2100,f)
# f.close()

# # # Output P45 2050
# with open(outdir + 'coben_P45_2050_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(coben_P45_2050_cell,f)
# f.close()

# with open(outdir + 'adapt_coben_P45_2050_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_coben_P45_2050_cell,f)
# f.close()

# with open(outdir + 'dif_coben_cell_P45_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(dif_coben_cell_P45_2050,f)
# f.close()

# # # Output P45 2100
# with open(outdir + 'coben_P45_2100_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(coben_P45_2100_cell,f)
# f.close()

# with open(outdir + 'adapt_coben_P45_2100_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_coben_P45_2100_cell,f)
# f.close()

# with open(outdir + 'dif_coben_cell_P45_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(dif_coben_cell_P45_2100,f)
# f.close()

# # Output P37 2050
# with open(outdir + 'coben_P37_2050_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(coben_P37_2050_cell,f)
# f.close()

# with open(outdir + 'adapt_coben_P37_2050_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_coben_P37_2050_cell,f)
# f.close()

# with open(outdir + 'dif_coben_cell_P37_2050_' + str(inf_rate), 'wb') as f:
#       cp.dump(dif_coben_cell_P37_2050,f)
# f.close()

# # # Output P37 2100
# with open(outdir + 'coben_P37_2100_cell_' + str(inf_rate), 'wb') as f: 
#       cp.dump(coben_P37_2100_cell,f)
# f.close()

# with open(outdir + 'adapt_coben_P37_2100_cell_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_coben_P37_2100_cell,f)
# f.close()

# with open(outdir + 'dif_coben_cell_P37_2100_' + str(inf_rate), 'wb') as f:
#       cp.dump(dif_coben_cell_P37_2100,f)
# f.close()

# # Output health burden data
# with open(outdir + 'base_burden_all_' + str(inf_rate), 'wb') as f:
#       cp.dump(base_burden,f)
# f.close()

# with open(outdir + 'adapt_burden_all_' + str(inf_rate), 'wb') as f:
#       cp.dump(adapt_burden,f)
# f.close()




