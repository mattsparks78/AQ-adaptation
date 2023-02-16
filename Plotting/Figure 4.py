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


# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir4 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Social Learning Model\\"
indir5 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Social Learning Model\Figure 4\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\Social Learning\\"

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


# Load pickle with net benefit amounts
with open(indir5 + 'net_ben_dict_sl_cost', 'rb') as f:
     nb = cp.load(f)
f.close()
    
# Select cases with lowest regret values
nb_opt = {}
# 2050 is 1.5, 0.5
nb_opt['2050'] = nb[1.5][0.5]
# 2100 is 0.9, 0.6
nb_opt['2100'] = nb[0.9][0.6]


# Load baseline case
with open(indir4 + 'net_ben_opt', 'rb') as f:
     nb_base = cp.load(f)
f.close()

# Calculate national mean values for best interventions
nb_opt_high = {}
nb_opt_low = {}
for i in nb_opt.keys():    
    nb_opt_high[i] = {}
    nb_opt_low[i] = {}
    for year in [2050, 2100]:
        nb_opt_high[i][year] = {}
        nb_opt_low[i][year] = {}
        for pol in pol_list:
            nat_list = []
            for IC in IC_list:
                scen_list = []
                for r, c in rc:
                    popu = pop.loc[(pop['ROW'] == r)& (pop['COL'] == c), str(year)].item()
                    scen_list.append(popu * np.mean(nb_opt[i][year][pol][IC][str((r, c))]))
                nat_list.append(sum(scen_list))
            nb_opt_high[i][year][pol] = max(nat_list) / 1E12
            nb_opt_low[i][year][pol] = min(nat_list) / 1E12
            
# Calculate national mean values for baseline adaptation case                
nb_base_high = {}
nb_base_low = {}
for year in [2050, 2100]:
    nb_base_high[year] = {}
    nb_base_low[year] = {}
    for pol in pol_list:
        nat_list = []
        for IC in IC_list:
            scen_list = []
            for r, c in rc:
                popu = pop.loc[(pop['ROW'] == r)& (pop['COL'] == c), str(year)].item()
                scen_list.append(popu * np.mean(nb_base[year][pol][IC][str((r, c))]))
            nat_list.append(sum(scen_list))
        nb_base_high[year][pol] = max(nat_list) / 1E12
        nb_base_low[year][pol] = min(nat_list) / 1E12

# Calculate national mean values for rational adaptation net benefits
rational_nb = {}
for pol in pol_list:
    rational_nb[pol] = {}
    fname = 'adapt_util_all_' + pol
    with open(indir1 + fname, 'rb') as f:
         ad_ut = cp.load(f)
    f.close()
    for year in [2050, 2100]:
        rational_nb[pol][year] = []
        if year == 2050:
            year_range = years[0]
        elif year == 2100:
            year_range = years[1]
        for Yr in year_range:
            for IC in IC_list:
                rational_nb[pol][year].append(sum(ad_ut[Yr][IC].values()))
        rational_nb[pol][year] = sum(rational_nb[pol][year]) / len(rational_nb[pol][year]) / 1E12

# Calculate series of 30 year net ben values for each IC
nb_30 = {}
for i in nb_opt.keys():
    nb_30[i] = {}
    for pol in pol_list:
        nb_30[i][pol] = {}
        for year in [2050, 2100]:
            nb_30[i][pol][year] = {}
            if year == 2050:
                year_range = range(2041, 2066)
            elif year == 2100:
                year_range = range(2091, 2116)
            for IC in IC_list:
                nb_30[i][pol][year][IC] = []
                for Yr in year_range:
                    nat_nb = []
                    for r, c in rc:
                        popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), str(year)].item()
                        key = str(Yr) + '-01-01 00:00:00'
                        nat_nb.append(popu * nb_opt[i][year][pol][IC][str((r, c))][key])
                    nb_30[i][pol][year][IC].append(sum(nat_nb) / 1E12)

# Calculate series of national values for rational nb
rat_nat_30 = {}
for pol in pol_list:
    rat_nat_30[pol] = {}
    fname = 'adapt_util_all_' + pol
    with open(indir1 + fname, 'rb') as f:
         ad_ut = cp.load(f)
    f.close()
    for year in [2050, 2100]:
        rat_nat_30[pol][year] = []
        if year == 2050:
            year_range = years[0]
        elif year == 2100:
            year_range = years[1]
        for Yr in year_range:
            for IC in IC_list:
                rat_nat_30[pol][year].append(sum(ad_ut[Yr][IC].values()) / 1E12)
                
# Calculate series of national values for no adaptation policy case
nb_30_base = {}
for pol in pol_list:
    nb_30_base[pol] = {}
    for year in [2050, 2100]:
        nb_30_base[pol][year] = {}
        if year == 2050:
            year_range = range(2041, 2066)
        elif year == 2100:
            year_range = range(2091, 2116)
        for IC in IC_list:
            nb_30_base[pol][year][IC] = []
            for Yr in year_range:
                nat_nb = []
                for r, c in rc:
                    popu = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), str(year)].item()
                    key = str(Yr) + '-01-01 00:00:00'
                    nat_nb.append(popu * nb_base[year][pol][IC][str((r, c))][key])
                nb_30_base[pol][year][IC].append(sum(nat_nb) / 1E12)

#%%
# Convert 2000 to 2020 $


#%%
### Plot

title = 'Social Learning Intervention Comparisons Bars Legend Only Opt Pol 2020'
# Plot all policy scenarios on the same plot PER CAPITA
fig, ax = plt.subplots(2, 3, figsize = (12, 10))

# Define general info
clist = ['black','olivedrab', 'mediumorchid', 'grey']
fsize = 16
lsize = 12
xlabels = ['High IC', 'Low IC']
width = 0.2
space = 0
msize2 = 48

# Plot cases
for i in [0, 1]:
    if i == 0:
        year = 2050
    elif i == 1:
        year = 2100
    for j in [0, 1, 2]:
        if j == 0:
            pol = 'REF'
        elif j == 1:
            pol = 'P45'
        elif j == 2:
            pol = 'P37'
    
        if i == 0:
            ax[i, j].set_xlim(0,3)
            
            ax[i, j].bar(1 - 0.5*width - space, round(nb_base_high[year][pol]*IF, 2), width = width, color = clist[0], label = 'No Adaptation Policy')
            ax[i, j].bar(1 + 0.5*width + space, round(nb_opt_high['2050'][year][pol]*IF, 2), width = width, color = clist[1], label = 'Adaptation Policy for 2050')
            
            
            ax[i, j].bar(2 - 0.5*width - space, round(nb_base_low[year][pol]*IF, 2),width = width, color = clist[0])
            ax[i, j].bar(2 + 0.5*width + space, round(nb_opt_low['2050'][year][pol]*IF, 2), width = width, color = clist[1])
            
        
            # First line and markers
            x = [1 - 0.5*width - space, 1 - 0.5*width - space]
            y = [np.min(nb_30_base[pol][year]['IC1'])*IF, np.max(nb_30_base[pol][year]['IC1'])*IF]
            ax[i, j].plot(x, y, color = clist[3])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            # Second line and markers
            x = [1 + 0.5*width + space, 1 + 0.5*width + space]
            y = [np.min(nb_30['2050'][pol][year]['IC1'])*IF, np.max(nb_30['2050'][pol][year]['IC1'])*IF]
            ax[i, j].plot(x, y, color = clist[0])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            
            # Fourth line and markers
            x = [2 - 0.5*width - space, 2 - 0.5*width - space]
            y = [np.min(nb_30_base[pol][year]['IC2'])*IF, np.max(nb_30_base[pol][year]['IC2'])*IF]
            ax[i, j].plot(x, y, color = clist[3])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            # Fifth line and markers
            x = [2 + 0.5*width + space, 2 + 0.5*width + space]
            y = [np.min(nb_30['2050'][pol][year]['IC2'])*IF, np.max(nb_30['2050'][pol][year]['IC2'])*IF]
            ax[i, j].plot(x, y, color = clist[0])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            
            
            
            
            ax[i, j].set_xticks([1, 2])
            
            sub_title = pol + ' ' + str(year)
            ax[i, j].set_title(sub_title, fontsize = fsize)
    
            ax[0, j].set_ylim(0, 12)
            ax[0, j].set_yticklabels([])
            ax[0, 0].set_yticklabels(list(range(0, 14,2)), fontsize = lsize)
            ax[0, j].set_xticklabels([])
            ax[1, j].set_ylim(0, 20)
            ax[1, j].set_yticklabels([])
            ax[1, 0].set_yticklabels(['0','2.5','5','7.5','10','12.5','15','17.5','20'], fontsize = lsize)
            ax[1, j].set_xticklabels(xlabels, fontsize = lsize)
            ax[i, 0].set_ylabel('Net Benefits of Adaptation\n$USD (2020) x 10^12', fontsize = fsize)
            ax[i, j].axhline(rational_nb[pol][year]*IF, color = 'grey', linestyle = '-', label = 'Optimal Adaptation Net Benefits')
            
        
        
        elif i == 1:
            
            ax[i, j].set_xlim(0,3)
            
            ax[i, j].bar(1 - 0.5*width - space, round(nb_base_high[year][pol]*IF, 2), width = width, color = clist[0], label = 'No Adaptation Policy')
            ax[i, j].bar(1 + 0.5*width + space, round(nb_opt_high['2100'][year][pol]*IF, 2),width = width, color = clist[2], label = 'Adaptation Policy for 2100')
            
            ax[i, j].bar(2 - 0.5*width - space, round(nb_base_low[year][pol]*IF, 2),width = width, color = clist[0])
            ax[i, j].bar(2 + 0.5*width + space,round(nb_opt_low['2100'][year][pol]*IF, 2), width = width, color = clist[2])
        
            # First line and markers
            x = [1 - 0.5*width - space, 1 - 0.5*width - space]
            y = [np.min(nb_30_base[pol][year]['IC2'])*IF, np.max(nb_30_base[pol][year]['IC2'])*IF]
            ax[i, j].plot(x, y, color = clist[3])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            
            # Third line and markers
            x = [1 + 0.5*width + space, 1 + 0.5*width + space]
            y = [np.min(nb_30['2100'][pol][year]['IC2'])*IF, np.max(nb_30['2100'][pol][year]['IC2'])*IF]
            ax[i, j].plot(x, y, color = clist[0])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[0])
            
            # Fourth line and markers
            x = [2 - 0.5*width - space, 2 - 0.5*width - space]
            y = [np.min(nb_30_base[pol][year]['IC3'])*IF, np.max(nb_30_base[pol][year]['IC3'])*IF]
            ax[i, j].plot(x, y, color = clist[3])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[3])
            
            # Sixth line and markers
            x = [2 + 0.5*width + space, 2 + 0.5*width + space]
            y = [np.min(nb_30['2100'][pol][year]['IC3'])*IF, np.max(nb_30['2100'][pol][year]['IC3'])*IF]
            ax[i, j].plot(x, y, color = clist[0])
            ax[i, j].scatter(x, y, marker = '_', s = msize2, color = clist[0])
            
            
            
            
            ax[i, j].set_xticks([1, 2])
            
            sub_title = pol + ' ' + str(year)
            ax[i, j].set_title(sub_title, fontsize = fsize)
    
            ax[0, j].set_ylim(0, 12)
            ax[0, j].set_xticklabels([])
            ax[1, j].set_ylim(0, 20)
            ax[1, j].set_xticklabels(xlabels, fontsize = lsize)
            ax[i, 0].set_ylabel('Net Benefits of Adaptation\n$USD (2020) x 10^12', fontsize = fsize)
            d = ax[i, j].axhline(rational_nb[pol][year]*IF, color = 'grey', linestyle = '-', label = 'Optimal Adaptation Net Benefits')
        
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

lines = [lines[0], lines[1], lines[2], lines[17]]
labels = [labels[0], labels[1], labels[2], labels[17]]
fig.legend(lines, labels, loc = 'lower center', fontsize = 16, ncol= 2)

plt.savefig(outdir + title + '.png', dpi = 400, bbox_inches = 'tight')
    
    

