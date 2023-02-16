# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:54:07 2022

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
import textwrap
from shapely.geometry import Polygon, MultiPolygon
import scipy as sp
import scipy.ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir4 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Social Learning Model\\"
indir5 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\US Boundary\\"
indir6 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\cb_2018_us_nation_5m\\'
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\\"

# Set general parameters

# Policy cases
pol_list = ['REF', 'P45', 'P37']

# IC cases
IC_list = ['IC'+ str(i) for i in range(1,6)]

# Population
pop = pd.read_csv(indir3 + 'pop_lepeule.csv')

# List of cells
rc = list(zip(pop['ROW'], pop['COL']))

# List of years
years = [range(2036,2066), range(2086, 2116)]

# Value of a Statistical Life
VSL = 7400000

# Economic Projection Factor
EPF = {}
for pol in pol_list:
    EPF[pol] = {}
EPF['REF'][2050] = 2.01
EPF['REF'][2100] = 3.94
EPF['P45'][2050] = 1.86
EPF['P45'][2100] = 3.50
EPF['P37'][2050] = 1.82
EPF['P37'][2100] = 3.37

# Population Projection Factor 
PPF = {}
PPF[2050] = 1.48
PPF[2100] = 1.73

#Conversion Factor from 2000 to 2020
IF = 1.4994

us_bound = gpd.read_file(indir6 + 'cb_2018_us_nation_5m.shp')
us_bound = us_bound.to_crs(epsg=4326)

# Calculate baseline national mortality timeseries
with open(indir1 + 'base_burden_all', 'rb') as f:
     base_val = cp.load(f)
f.close()

base_mort = {}

for pol in pol_list:
    base_mort[pol] = {}
    for year_set in years:
        for Yr in year_set:
            base_mort[pol][Yr] = {}
            if Yr < 2070:
                pol_year = 2050
            else:
                pol_year = 2100
            for IC in IC_list:
                base_mort[pol][Yr][IC] = sum(base_val[pol][Yr][IC].values()) / (VSL * EPF[pol][pol_year]) * PPF[pol_year]

# Calculate social learning national mortality timeseries
sl_mort = {}
for pol in pol_list:
    sl_mort[pol] = {}
    for pol_year in [2050, 2100]:
        sl_mort[pol][pol_year] = {}
        for IC in IC_list:
            sl_mort[pol][pol_year][IC] = {}
            file = 'mort_pc_' + pol + '_' + str(pol_year) + '_' + IC + '_density.csv'
            mort_pc = pd.read_csv(indir4 + file)
            mort_pc = mort_pc.to_dict()
            for i in range(0, 100, 4):
                mort_tot = []
                if pol_year == 2050:
                    Yr = int(2041 + i/4)
                else:
                    Yr = int(2091 + i/4)
                for r, c in rc:
                    percap = mort_pc[str((r, c))][i]
                    popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), str(pol_year)].item()
                    mort_tot.append(percap * popu)
                    
                sl_mort[pol][pol_year][IC][Yr] = sum(mort_tot)
                

# Calculate optimal national mortality timeseries
with open(indir1 + 'adapt_burden_all', 'rb') as f:
     opt_val = cp.load(f)
f.close()

opt_mort = {}

for pol in pol_list:
    opt_mort[pol] = {}
    for year_set in years:
        for Yr in year_set:
            opt_mort[pol][Yr] = {}
            if Yr < 2070:
                pol_year = 2050
            else:
                pol_year = 2100
            for IC in IC_list:
                opt_mort[pol][Yr][IC] = sum(opt_val[pol][Yr][IC].values()) / (VSL * EPF[pol][pol_year]) * PPF[pol_year]



# Calculate social learning national cost timeseries
sl_cost = {}
for pol in pol_list: 
    sl_cost[pol] = {}
    for pol_year in [2050, 2100]:
        sl_cost[pol][pol_year] = {}
        for IC in IC_list:
            sl_cost[pol][pol_year][IC] = {}
            file = 'cost_pc_' + pol + '_' + str(pol_year) + '_' + IC + '_density.csv'
            cost_pc = pd.read_csv(indir4 + file)
            cost_pc = cost_pc.to_dict()
            for i in range(0, 100, 4):
                cost_tot = []
                if pol_year == 2050:
                    Yr = int(2041 + i/4)
                else:
                    Yr = int(2091 + i/4)
                for r, c in rc:
                    percap = cost_pc[str((r, c))][i]
                    popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), str(pol_year)].item()
                    cost_tot.append(percap * popu)
                    
                sl_cost[pol][pol_year][IC][Yr] = sum(cost_tot)


# Calculate optimal national cost timeseries
opt_cost = {}

for pol in pol_list:
    opt_cost[pol] = {}
    file = 'adapt_cost_all_' + pol
    with open(indir1 + file, 'rb') as f:
         opt_cost_pol = cp.load(f)
    f.close()
    for year_set in years:
        for Yr in year_set:
            opt_cost[pol][Yr] = {}
            for IC in IC_list:
                opt_cost[pol][Yr][IC] = sum(opt_cost_pol[Yr][IC].values())


# Plot all policy scenarios on the same plot 
# for pol_year in [2050, 2100]:
#     if pol_year == 2050:
#         year_range = range(2041, 2066)
#     else:
#         year_range = range(2091, 2116)
#     base_mort_list = {}
#     sl_mort_list = {}
#     opt_mort_list = {}
#     sl_cost_list = {}
#     opt_cost_list = {}
#     for pol in pol_list:
#         base_mort_list[pol] = []
#         sl_mort_list[pol] = []
#         opt_mort_list[pol] = []
#         sl_cost_list[pol] = []
#         opt_cost_list[pol] = []
#         for Yr in year_range:
#             base_mort_list[pol].append(sum(base_mort[pol][Yr].values()) / 5)
#             sl_IC_mort_list = []
#             opt_mort_list[pol].append(sum(opt_mort[pol][Yr].values()) / 5)
#             sl_IC_cost_list = []
#             opt_cost_list[pol].append((sum(opt_cost[pol][Yr].values()) / 5) / 1E12)
#             for IC in IC_list:
#                 sl_IC_mort_list.append(sl_mort[pol][pol_year][IC][Yr])
#                 sl_IC_cost_list.append(sl_cost[pol][pol_year][IC][Yr])
#             sl_mort_list[pol].append(np.mean(sl_IC_mort_list))
#             sl_cost_list[pol].append(np.mean(sl_IC_cost_list) / 1E12)
    
#     title = 'Social Learning vs Optimal ' + str(pol_year)

#     fig, axs = plt.subplots(2)
#     # ax.plot(list(year_range), base_mort_list, color = 'blue', marker = 'o')
#     line1, = axs[0].plot(list(year_range), opt_mort_list['P37'], color = 'skyblue', linestyle = '-.', linewidth = 0.75)
#     line2, = axs[0].plot(list(year_range), opt_mort_list['P45'], color = 'blue',  linestyle = '--', linewidth = 0.75)
#     line3, = axs[0].plot(list(year_range), opt_mort_list['REF'], color = 'navy', linewidth = 0.75)
      
#     line4, = axs[0].plot(list(year_range), sl_mort_list['P37'], color = 'darksalmon', linestyle = '-.', linewidth = 0.75)
#     line5, = axs[0].plot(list(year_range), sl_mort_list['P45'], color = 'red', linestyle = '--', linewidth = 0.75)
#     line6, = axs[0].plot(list(year_range), sl_mort_list['REF'], color = 'maroon', linewidth = 0.75)
    
#     axs[1].plot(list(year_range), opt_cost_list['P37'], color = 'skyblue', linestyle = '-.', linewidth = 0.75)
#     axs[1].plot(list(year_range), opt_cost_list['P45'], color = 'blue', linestyle = '--', linewidth = 0.75)
#     axs[1].plot(list(year_range), opt_cost_list['REF'], color = 'navy', linewidth = 0.75)
    
#     axs[1].plot(list(year_range), sl_cost_list['P37'], color = 'darksalmon', linestyle = '-.', linewidth = 0.75)
#     axs[1].plot(list(year_range), sl_cost_list['P45'], color = 'red', linestyle = '--', linewidth = 0.75)
#     axs[1].plot(list(year_range), sl_cost_list['REF'], color = 'maroon', linewidth = 0.75)
    
#     axs[1].set_xlabel('Year')
#     axs[0].set_ylabel('Mortalities')
#     axs[1].set_ylabel('Adaptation Cost ($USD2000x10^12)')
#     axs[1].yaxis.set_label_coords(-0.135, 0.4)
#     axs[0].set_title(title)
#     mpl.legend_handler.HandlerTuple(ndivide = 2)
#     axs[0].legend([line6, line3],['Social Learning', 'Optimal'])
#     axs[1].legend([(line3, line6), (line2, line5), (line1, line4)],['REF', 'P45', 'P37'])
#     axs[0].axes.get_xaxis().set_ticklabels([])
#     axs[0].set_ylim([0, None])
#     axs[1].set_ylim([0, None])
    
    
  
for pol_year in [2050]:
    if pol_year == 2050:
        year_range = range(2041, 2066)
        popu = 267818227
    else:
        year_range = range(2091, 2116)
        popu = 313057802
    base_mort_list = {}
    sl_mort_list = {}
    opt_mort_list = {}
    sl_cost_list = {}
    opt_cost_list = {}
    for pol in pol_list:
        base_mort_list[pol] = []
        sl_mort_list[pol] = []
        opt_mort_list[pol] = []
        sl_cost_list[pol] = []
        opt_cost_list[pol] = []
        for Yr in year_range:
            base_mort_list[pol].append((sum(base_mort[pol][Yr].values()) / 5) / popu)
            sl_IC_mort_list = []
            opt_mort_list[pol].append((sum(opt_mort[pol][Yr].values()) / 5) / popu)
            sl_IC_cost_list = []
            opt_cost_list[pol].append((sum(opt_cost[pol][Yr].values()) / 5) / popu)
            for IC in IC_list:
                sl_IC_mort_list.append(sl_mort[pol][pol_year][IC][Yr])
                sl_IC_cost_list.append(sl_cost[pol][pol_year][IC][Yr])
            sl_mort_list[pol].append(np.mean(sl_IC_mort_list) / popu)
            sl_cost_list[pol].append(np.mean(sl_IC_cost_list) / popu)
          

# Add bar chart of adaptation cost per avoided mortality
# Calculate total mortalities under baseline adaptation in 2050
with open(indir1 + 'base_burden_all', 'rb') as f:
     bba = cp.load(f)
f.close()

# Set economic factors for 2050

EF_2050 = {'REF':2.01 , 'P45':1.86, 'P37':1.82 }

vals_b = {}
morts_b = {}
for pol in ['REF', 'P45', 'P37']:
    vals_b[pol] = []
    for year in range(2036, 2066):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            vals_b[pol].append(sum(bba[pol][year][IC].values()))

    mean_vals = sum(vals_b[pol]) / len(vals_b[pol])
    # Calculate mortalities by valuation total/ valuation per. Multiplying VSL by economic factor for REF 2050
    morts_b[pol] = mean_vals / (7400000 * EF_2050[pol])

# Calculate total mortalities under optimal adaptation in 2050
with open(indir1 + 'adapt_burden_all', 'rb') as f:
     aba = cp.load(f)
f.close()

vals_ad = {}
morts_ad = {}
for pol in ['REF', 'P45', 'P37']:
    vals_ad[pol] = []
    for year in range(2036, 2066):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            vals_ad[pol].append(sum(aba[pol][year][IC].values()))

    mean_vals = sum(vals_ad[pol]) / len(vals_ad[pol])
    # Calculate mortalities by valuation total/ valuation per. Multiplying VSL by economic factor for REF 2050
    morts_ad[pol] = mean_vals / (7400000 * EF_2050[pol])

# Calculate total cost of optimal adaptation in 2050
costs_opt = {}
for pol in ['REF', 'P45', 'P37']:
    costs_opt[pol] = []
    fname = 'adapt_cost_all_' + pol
    with open(indir1 + fname, 'rb') as f:
         cost_dict = cp.load(f)
    f.close()
    for year in range(2036, 2066):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            costs_opt[pol].append(sum(cost_dict[year][IC].values()))
    costs_opt[pol] = sum(costs_opt[pol]) / len(costs_opt[pol])
    
# Calculate total mortalities under social learning in 2050
morts_sl = {}
for pol in ['REF', 'P45', 'P37']:
    morts_sl[pol] = []
    for IC in ['IC' + str(i) for i in range(1,6)]:
        fname = 'mort_pc_' + pol + '_2050_' + IC + '_density.csv' 
        mort_df = pd.read_csv(indir4 + fname)    
        for i in range(0, 99, 4):
            a = mort_df.iloc[i]
            mort_list = []
            for r, c in rc:
                mort_rate = a[str((r, c))]
                popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2050'].item()
                tot_mort = popu * mort_rate
                mort_list.append(tot_mort)
            morts_sl[pol].append(sum(mort_list))
    morts_sl[pol] = sum(morts_sl[pol]) / len(morts_sl[pol])

# Calculate total cost under social learning in 2050
costs_sl = {}
for pol in ['REF', 'P45', 'P37']:
    costs_sl[pol] = []
    for IC in ['IC' + str(i) for i in range(1,6)]:
        fname = 'cost_pc_' + pol + '_2050_' + IC + '_density.csv' 
        cost_df = pd.read_csv(indir4 + fname)    
        for i in range(0, 99, 4):
            a = cost_df.iloc[i]
            cost_list = []
            for r, c in rc:
                cost_rate = a[str((r, c))]
                popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2050'].item()
                tot_cost = popu * cost_rate
                cost_list.append(tot_cost)
            costs_sl[pol].append(sum(cost_list))
    costs_sl[pol] = sum(costs_sl[pol]) / len(costs_sl[pol])

# Calculate cost effectiveness of optimal and social learning in 2050
ce = {}
ce['SL'] = {}
ce['Opt'] = {}
for pol in ['REF', 'P45', 'P37']:
    ce['SL'][pol] = costs_sl[pol] / (morts_b[pol] - morts_sl[pol]) 
    ce['Opt'][pol] = costs_opt[pol] / (morts_b[pol] - morts_ad[pol]) 


# Calculate total mortalities under baseline adaptation in 2100

EF_2100 = {'REF':3.94 , 'P45':3.5, 'P37':3.37 }

vals_b_2100 = {}
morts_b_2100 = {}
for pol in ['REF', 'P45', 'P37']:
    vals_b_2100[pol] = []
    for year in range(2086, 2116):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            vals_b_2100[pol].append(sum(bba[pol][year][IC].values()))

    mean_vals = sum(vals_b_2100[pol]) / len(vals_b_2100[pol])
    # Calculate mortalities by valuation total/ valuation per. Multiplying VSL by economic factor for REF 2050
    morts_b_2100[pol] = mean_vals / (7400000 * EF_2100[pol])

# Calculate total mortalities under optimal adaptation in 2050

vals_ad_2100 = {}
morts_ad_2100 = {}
for pol in ['REF', 'P45', 'P37']:
    vals_ad_2100[pol] = []
    for year in range(2086, 2116):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            vals_ad_2100[pol].append(sum(aba[pol][year][IC].values()))

    mean_vals = sum(vals_ad[pol]) / len(vals_ad[pol])
    # Calculate mortalities by valuation total/ valuation per. Multiplying VSL by economic factor for REF 2050
    morts_ad_2100[pol] = mean_vals / (7400000 * EF_2100[pol])

# Calculate total cost of optimal adaptation in 2050
costs_opt_2100 = {}
for pol in ['REF', 'P45', 'P37']:
    costs_opt_2100[pol] = []
    fname = 'adapt_cost_all_' + pol
    with open(indir1 + fname, 'rb') as f:
         cost_dict = cp.load(f)
    f.close()
    for year in range(2086, 2116):
        for IC in ['IC' + str(i) for i in range(1,6)]:
            costs_opt_2100[pol].append(sum(cost_dict[year][IC].values()))
    costs_opt_2100[pol] = sum(costs_opt_2100[pol]) / len(costs_opt_2100[pol])
    
# Calculate total mortalities under social learning in 2050
morts_sl_2100 = {}
for pol in ['REF', 'P45', 'P37']:
    morts_sl_2100[pol] = []
    for IC in ['IC' + str(i) for i in range(1,6)]:
        fname = 'mort_pc_' + pol + '_2100_' + IC + '_density.csv' 
        mort_df = pd.read_csv(indir4 + fname)    
        for i in range(0, 99, 4):
            a = mort_df.iloc[i]
            mort_list = []
            for r, c in rc:
                mort_rate = a[str((r, c))]
                popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2100'].item()
                tot_mort = popu * mort_rate
                mort_list.append(tot_mort)
            morts_sl_2100[pol].append(sum(mort_list))
    morts_sl_2100[pol] = sum(morts_sl_2100[pol]) / len(morts_sl_2100[pol])

# Calculate total cost under social learning in 2050
costs_sl_2100 = {}
for pol in ['REF', 'P45', 'P37']:
    costs_sl_2100[pol] = []
    for IC in ['IC' + str(i) for i in range(1,6)]:
        fname = 'cost_pc_' + pol + '_2100_' + IC + '_density.csv' 
        cost_df = pd.read_csv(indir4 + fname)    
        for i in range(0, 99, 4):
            a = cost_df.iloc[i]
            cost_list = []
            for r, c in rc:
                cost_rate = a[str((r, c))]
                popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2050'].item()
                tot_cost = popu * cost_rate
                cost_list.append(tot_cost)
            costs_sl_2100[pol].append(sum(cost_list))
    costs_sl_2100[pol] = sum(costs_sl_2100[pol]) / len(costs_sl_2100[pol])

# Calculate cost effectiveness of optimal and social learning in 2050
ce_2100 = {}
ce_2100['SL'] = {}
ce_2100['Opt'] = {}
for pol in ['REF', 'P45', 'P37']:
    ce_2100['SL'][pol] = costs_sl_2100[pol] / (morts_b_2100[pol] - morts_sl_2100[pol]) 
    ce_2100['Opt'][pol] = costs_opt_2100[pol] / (morts_b_2100[pol] - morts_ad_2100[pol]) 

#%%
#Calculate regional proportion of adapters timeseries
# Load grid cell contribution to region file
grid_reg = pd.read_csv(indir2 + 'grid_by_census_region.csv')

# Load region population file
with open(indir2 + 'region_pop_lepeule', 'rb') as f:
     pop_reg = cp.load(f)
f.close()

# Calculate proportion of regional population in each grid cell
reg_prop = grid_reg[['ROW', 'COL']]
for reg in ['Northeast', 'Midwest', 'South', 'West']:
    reg_prop[reg] = 0
    s = grid_reg[grid_reg[reg] != 0]
    rc_temp = list(zip(s['ROW'], s['COL']))
    for r, c in rc_temp:
        prop = s.loc[(s['ROW'] == r) & (s['COL'] == c), reg].item()
        popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item()
        pop_cell_reg = popu * prop
        pop_prop = pop_cell_reg / pop_reg[reg]
        reg_prop.loc[(reg_prop['ROW'] == r) & (reg_prop['COL'] == c), reg] = pop_prop

# # Test
# for reg in ['Northeast', 'Midwest', 'South', 'West']:
#     print(sum(reg_prop[reg]))

# Calculate xi time series by region for REF 2050
xi_REF_2050 = {}
for IC in IC_list:
    xi_REF_2050[IC] = {}
    file = 'x_' + 'REF' + '_' + str(2050) + '_' + IC + '_density.csv'
    xi_vals = pd.read_csv(indir4 + file)
    xi_vals = xi_vals.to_dict()
    for reg in ['Northeast', 'Midwest', 'South', 'West']:
        xi_REF_2050[IC][reg] = []
        for i in range(0, 100, 4):
            xi_reg = []
            for r, c in rc:
                xi_cell = xi_vals[str((r, c))][i]
                portion_reg = reg_prop.loc[(reg_prop['ROW'] == r) & (reg_prop['COL'] == c), reg].item()
                xi_reg.append(xi_cell * portion_reg)
            xi_REF_2050[IC][reg].append(sum(xi_reg))     
        
xi_ts = {}
for reg in ['Northeast', 'Midwest', 'South', 'West']:
    xi_ts[reg] = []
    for i in range(len(xi_REF_2050[IC][reg])):
        reg_vals = []
        for IC in IC_list:
            reg_vals.append(xi_REF_2050[IC][reg][i])
        xi_ts[reg].append(sum(reg_vals) / 5)

#%%
#Create geodataframes for adaptation utility by region for RA and SL models
# Rational Actor Model REF
RA_util = gpd.read_file(indir2 + 'US_census_regions.shp')
RA_util['Utility'] = 0

# Calculate regional adaptation utility for REF 2050
# Load adaptation utility file
with open(indir1 + 'adapt_util_pc_all_REF', 'rb') as f:
     RA_ad_util = cp.load(f)
f.close()

# Adaptation utility by cell
RA_cell = {}
for r, c in rc:
    RA_cell[(r, c)] = []
    for year in range(2036,2066):
        for IC in IC_list:
            RA_cell[(r, c)].append(RA_ad_util[year][IC][(r, c)])
    RA_cell[(r, c)] = np.mean(RA_cell[(r, c)])
    
# Adaptation utility by region
for reg in ['Northeast', 'Midwest', 'South', 'West']:
    reg_util = []
    for r, c in rc:
        ad_util_cell = RA_cell[(r, c)]
        portion_reg = reg_prop.loc[(reg_prop['ROW'] == r) & (reg_prop['COL'] == c), reg].item()
        reg_util.append(ad_util_cell * portion_reg)
    RA_util.loc[RA_util['regions'] == reg, 'Utility'] = sum(reg_util)
    
# Rational Actor Model P37
RA_util_P37 = gpd.read_file(indir2 + 'US_census_regions.shp')
RA_util_P37['Utility'] = 0

# Calculate regional adaptation utility for REF 2050
# Load adaptation utility file
with open(indir1 + 'adapt_util_pc_all_P37', 'rb') as f:
     RA_ad_util_P37 = cp.load(f)
f.close()

# Adaptation utility by cell
RA_cell_P37 = {}
for r, c in rc:
    RA_cell_P37[(r, c)] = []
    for year in range(2036,2066):
        for IC in IC_list:
            RA_cell_P37[(r, c)].append(RA_ad_util_P37[year][IC][(r, c)])
    RA_cell_P37[(r, c)] = np.mean(RA_cell_P37[(r, c)])
    
RA_diff = {}
for r, c in rc:
    RA_diff[(r, c)] = RA_cell[(r, c)] - RA_cell_P37[(r, c)] 


#%%
# Social Learning Model
SL_util_cell = {}

for r, c in rc:
    SL_util_cell[(r, c)] = []
    for IC in IC_list:
        SL_util = pd.read_csv(indir4 + 'net_ben_REF_2050_' + IC +'_density.csv')
        SL_util = SL_util.drop(columns = ['Unnamed: 0'])
        SL_util_cell[(r, c)].append(SL_util[str((r, c))].mean())
    SL_util_cell[(r, c)] = sum(SL_util_cell[(r, c)]) / len(SL_util_cell[(r, c)])
    
# Social Learning Model P37
SL_util_cell_P37 = {}

for r, c in rc:
    SL_util_cell_P37[(r, c)] = []
    for IC in IC_list:
        SL_util_P37 = pd.read_csv(indir4 + 'net_ben_P37_2050_' + IC +'_density.csv')
        SL_util_cell_P37[(r, c)].append(SL_util_P37[str((r, c))].mean())
    SL_util_cell_P37[(r, c)] = sum(SL_util_cell_P37[(r, c)]) / len(SL_util_cell_P37[(r, c)])
    
# Calc diff
SL_diff = {}
for r, c in rc:
    SL_diff[(r, c)] = SL_util_cell[(r, c)] - SL_util_cell_P37[(r, c)] 

    
#%% Create contour plots for RA and SL

# Plot Rational Actor Map
# Import grid  gdf
grid = gpd.read_file(indir2 + "BenMapGridPoints.csv")
grid = grid.drop('geometry', axis = 1)
grid['geom'] = grid['geom'].apply(wkt.loads)
grid = grid.set_geometry('geom')
grid = grid.set_crs(epsg=4326)
grid['ROW'] = grid['ROW'].astype('int')
grid['COL'] = grid['COL'].astype('int')
grid = grid.drop([1, 219, 220])

# Plot adaptation days by grid cell
RA_gdf = copy.deepcopy(grid)
RA_gdf['pc_utility'] = 0

for r, c in rc:
    RA_gdf.loc[(RA_gdf['ROW'] == r) & (RA_gdf['COL'] == c), 'pc_utility'] = RA_diff[(r, c)]
    
SL_gdf = copy.deepcopy(grid)
SL_gdf['pc_utility'] = 0

for r, c in rc:
    SL_gdf.loc[(SL_gdf['ROW'] == r) & (SL_gdf['COL'] == c), 'pc_utility'] = SL_diff[(r, c)]
    
#Convert pc utility to 2020
RA_gdf['pc_utility'] = RA_gdf['pc_utility'] * IF
SL_gdf['pc_utility'] = SL_gdf['pc_utility'] * IF
    

mymin = -12000
mymax = 12000

norm = mpl.colors.TwoSlopeNorm(vmin= mymin, vcenter = 0, vmax= mymax)

# Plot Contours
# Define function to convert contour to gdf
def collec_to_gdf(collec_poly):
     #"""Transform a `matplotlib.contour.QuadContourSet` to a GeoDataFrame"""
     polygons, colors = [], []
     for i, polygon in enumerate(collec_poly.collections):
         mpoly = []
         for path in polygon.get_paths():
             try:
                 path.should_simplify = False
                 poly = path.to_polygons()
                 # Each polygon should contain an exterior ring + maybe hole(s):
                 exterior, holes = [], []
                 if len(poly) > 0 and len(poly[0]) > 3:
                     # The first of the list is the exterior ring :
                     exterior = poly[0]
                     # Other(s) are hole(s):
                     if len(poly) > 1:
                         holes = [h for h in poly[1:] if len(h) > 3]
                 mpoly.append(Polygon(exterior, holes).buffer(0))
             except:
                 print('Warning: Geometry error when making polygon #{}'
                       .format(i))
         if len(mpoly) > 1:
             mpoly = MultiPolygon(mpoly)
             polygons.append(mpoly)
             colors.append(polygon.get_facecolor().tolist()[0])
         elif len(mpoly) == 1:
             polygons.append(mpoly[0])
             colors.append(polygon.get_facecolor().tolist()[0])
     return gpd.GeoDataFrame(
         geometry=polygons,
         data={'RGBA': colors},
         crs = 'epsg:4326')


# Create array of x and y values
lon_arr = grid['Lon'].unique()
lat_arr = grid['Lat'].unique()

# Change datatype of lon/lat arrays
lon_arr = lon_arr.astype('float')
lat_arr = lat_arr.astype('float')

#Sort arrays
lon_arr = np.sort(lon_arr)
lat_arr = np.sort(lat_arr)

# Create empty array to store pol ben values
ben_arr = np.zeros((len(lat_arr), len(lon_arr)))

for i in range(len(lat_arr)):
    for j in range(len(lon_arr)):
        b = RA_gdf.loc[(RA_gdf['Lon'] == str(lon_arr[j])) & (RA_gdf['Lat'] == str(lat_arr[i])), 'pc_utility']
        if len(b) == 0:
            pass 
        else:
            ben_arr[i, j] = RA_gdf.loc[(RA_gdf['Lon'] == str(lon_arr[j])) & (RA_gdf['Lat'] == str(lat_arr[i])), 'pc_utility'].item() 

lon_arr = sp.ndimage.zoom(lon_arr, 5)
lat_arr = sp.ndimage.zoom(lat_arr, 5)
ben_arr = sp.ndimage.zoom(ben_arr, 5)

        
collec_poly = plt.contourf(lon_arr, lat_arr, ben_arr, levels = 25, cmap = 'coolwarm', norm = norm);

test = collec_to_gdf(collec_poly)

ticks=[mymin, mymin*2/3,mymin/3,0,mymax/3,mymax*2/3,mymax]
label = 'Adaptation Utility ($USD 2020)'
tick_labels = ['-$'+str(mymax), '-$'+str(mymax*2/3),'-$'+str(mymax/3),'$0','$'+str(mymax/3),'$'+str(mymax*2/3),'$'+str(mymax)]


colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
test['RGBA'] = colors
test = test.to_crs(epsg=4326)
for i in range(len(test.geometry)):
    test['geometry'][i] = test['geometry'][i].intersection(us_bound['geometry'][0])


###### Social Learning Model Plot ############

# Create array of x and y values
lon_arr = grid['Lon'].unique()
lat_arr = grid['Lat'].unique()

# Change datatype of lon/lat arrays
lon_arr = lon_arr.astype('float')
lat_arr = lat_arr.astype('float')

#Sort arrays
lon_arr = np.sort(lon_arr)
lat_arr = np.sort(lat_arr)

# Create empty array to store pol ben values
ben_arr = np.zeros((len(lat_arr), len(lon_arr)))

for i in range(len(lat_arr)):
    for j in range(len(lon_arr)):
        b = SL_gdf.loc[(SL_gdf['Lon'] == str(lon_arr[j])) & (SL_gdf['Lat'] == str(lat_arr[i])), 'pc_utility']
        if len(b) == 0:
            pass 
        else:
            ben_arr[i, j] = SL_gdf.loc[(SL_gdf['Lon'] == str(lon_arr[j])) & (SL_gdf['Lat'] == str(lat_arr[i])), 'pc_utility'].item() 

lon_arr = sp.ndimage.zoom(lon_arr, 5)
lat_arr = sp.ndimage.zoom(lat_arr, 5)
ben_arr = sp.ndimage.zoom(ben_arr, 5)

        
collec_poly = plt.contourf(lon_arr, lat_arr, ben_arr, levels = 25, cmap = 'coolwarm', norm = norm);


SL_plot = collec_to_gdf(collec_poly)

ticks=[mymin, mymin*2/3,mymin/3,0,mymax/3,mymax*2/3,mymax]
label = 'Net benefits of adaptation per capita\nin 2050 under reference vs. climate policy\n(REF vs P37) ($USD 2020)'
tick_labels = ['-$'+str(round(mymax)), '-$'+str(round(mymax*2/3)),'-$'+str(round(mymax/3)),'$0','$'+str(round(mymax/3)),'$'+str(round(mymax*2/3)),'$'+str(round(mymax))]


colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
SL_plot['RGBA'] = colors
SL_plot = SL_plot.to_crs(epsg=4326)
for i in range(len(SL_plot.geometry)):
    SL_plot['geometry'][i] = SL_plot['geometry'][i].intersection(us_bound['geometry'][0])


# Convert SL and RA costs to 2020
for pol in pol_list:
    for i in range(0,25):
        sl_cost_list[pol][i] = sl_cost_list[pol][i] * IF
        opt_cost_list[pol][i] = opt_cost_list[pol][i] * IF


#%%
### Plot


title = 'SL vs RA Per Capita TS and Regional'
# Plot all policy scenarios on the same plot PER CAPITA
fig = plt.figure(figsize=(10, 8))

subfigs = fig.subfigures(1, 2, wspace = 0,hspace = 0.25, width_ratios = [1,1.4])

axsLeft = subfigs[0].subplots(3, 1)
axsRight = subfigs[1].subplots(3, 1, gridspec_kw = {'height_ratios':[1,1,0.12]})


# plt.subplots_adjust(wspace= 0.25, hspace= 0.25)

ax1 = axsLeft[0]
ax2 = axsLeft[1]
ax3 = axsLeft[2]
ax4 = axsRight[0]
ax5 = axsRight[1]
ax6 = axsRight[2]

# ax.plot(list(year_range), base_mort_list, color = 'blue', marker = 'o')
line1, = ax1.plot(list(year_range), opt_mort_list['P37'], color = 'skyblue', linestyle = '-.', linewidth = 0.75)
line2, = ax1.plot(list(year_range), opt_mort_list['P45'], color = 'blue',  linestyle = '--', linewidth = 0.75)
line3, = ax1.plot(list(year_range), opt_mort_list['REF'], color = 'navy', linewidth = 0.75)
  
line4, = ax1.plot(list(year_range), sl_mort_list['P37'], color = 'darksalmon', linestyle = '-.', linewidth = 0.75)
line5, = ax1.plot(list(year_range), sl_mort_list['P45'], color = 'red', linestyle = '--', linewidth = 0.75)
line6, = ax1.plot(list(year_range), sl_mort_list['REF'], color = 'maroon', linewidth = 0.75)

ax2.plot(list(year_range), opt_cost_list['P37'], color = 'skyblue', linestyle = '-.', linewidth = 0.75)
ax2.plot(list(year_range), opt_cost_list['P45'], color = 'blue', linestyle = '--', linewidth = 0.75)
ax2.plot(list(year_range), opt_cost_list['REF'], color = 'navy', linewidth = 0.75)

ax2.plot(list(year_range), sl_cost_list['P37'], color = 'darksalmon', linestyle = '-.', linewidth = 0.75)
ax2.plot(list(year_range), sl_cost_list['P45'], color = 'red', linestyle = '--', linewidth = 0.75)
ax2.plot(list(year_range), sl_cost_list['REF'], color = 'maroon', linewidth = 0.75)

line7, = ax3.plot(list(year_range), xi_ts['Northeast'], color = 'indigo', linestyle = '-.', linewidth = 0.75)
line8, = ax3.plot(list(year_range), xi_ts['Midwest'], color = 'rebeccapurple', linestyle = '--', linewidth = 0.75)
line9, = ax3.plot(list(year_range), xi_ts['South'], color = 'mediumorchid', linestyle = 'dotted', linewidth = 0.75)
line10, = ax3.plot(list(year_range), xi_ts['West'], color = 'fuchsia', linewidth = 0.75)


ylabel1 = "Mortality Rate\nPer Capita"
ylabel2 = 'Adaptation Cost\nPer Capita\n$USD(2020)'
ylabel3 = 'Adapter Proportion\nBy Census Region\nREF 2050'

fsize = 14

ax3.set_xlabel('Year', fontsize = fsize)
ax2.set_ylabel(ylabel2, fontsize = fsize)
ax1.set_ylabel(ylabel1, fontsize = fsize)
mpl.legend_handler.HandlerTuple(ndivide = 2)
ax1.legend([line6, line3],['Socially Influenced', 'Optimal'], fontsize = 12)
ax2.legend([(line3, line6), (line2, line5), (line1, line4)],['REF', 'P45', 'P37'], fontsize = 12, ncol = 3)

ax1.set_ylim([0, 0.002])
ax2.set_ylim([0, 25000])
ax2.set_xlim(2041, 2066)
ax1.annotate('(a)', (2063, 0.00015), fontsize = 14)
ax2.annotate('(b)', (2063, 22000), fontsize = 14)
ax3.annotate('(c)', (2063, 0.15), fontsize = 14)
ax4.annotate('(d)',(-120,27), fontsize = 14)
ax5.annotate('(e)',(-120,27), fontsize = 14)

ax3.set_ylim([0,1])
ax3.set_xlim(2041,2066)


ax3.set_ylabel(ylabel3, fontsize = fsize)
ax3.legend([line7, line8, line9, line10], ['Northeast','Midwest','South','West'], fontsize = 12, ncol = 2, loc = 'upper left')
ax3.annotate('(c)', (0.2, 2), fontsize = 14)

test.plot(color = test['RGBA'], ax = ax4)



# Set axes off
ax4.set_axis_off()
ax4.set_title('Optimal', fontsize = fsize)
# Plot regional boundaries
grid_reg = gpd.read_file(indir2 + 'US_census_regions.shp')

grid_reg.boundary.plot(ax = ax4, color = 'black')



SL_plot.plot(color = SL_plot['RGBA'], ax = ax5)


# Set axes off
ax5.set_axis_off()
ax5.set_title('Socially Influenced', fontsize = fsize)
grid_reg.boundary.plot(ax = ax5, color = 'black')

# Plot colorbar

cbar = mpl.colorbar.ColorbarBase(cmap = 'coolwarm', orientation='horizontal', ticks = ticks, norm = norm, ax = ax6)
cbar.set_label(label, fontsize = fsize)
cbar.ax.set_xticklabels(tick_labels)
# ax6.set_axis_off()


plt.savefig(outdir + '5_panel_fig_3_contour_2020.png', dpi = 400, bbox_inches = 'tight')