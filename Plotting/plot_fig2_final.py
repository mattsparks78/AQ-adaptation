# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 14:35:04 2023

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
from shapely.geometry import Polygon, MultiPolygon
import scipy as sp
import scipy.ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

# Set directories
indir1 = r".\Data\Sensitivity\\"
indir2 = r".\Data\Historical AQI\\"
indir3 = r".\Data\General\\"
indir4 = r".\Data\Social Learning Model\\"
indir5 = r".\Data\Historical Model PM\\"
indir6 = r'.\Data\cb_2018_us_nation_5m\\'
indir7 = r'.\Data\BenMap Grid\\'
indir8 = r".\Data\Future AQ\\"
indir9 = r".\Data\Valuations All Infiltration\\"
indir10 = r".\Data\Threshold\\"
outdir = r".\Figures\\"


# Population
pop = pd.read_csv(indir3 + 'pop_lepeule_race.csv')

# List of cells
rc = list(zip(pop['ROW'], pop['COL']))

#%%
# Load net benefits files
# Climate variability
with open(indir1 + 'climate_variation_net_benefits_for_conc_1hr', 'rb') as f:
     climate_nb = cp.load(f)
f.close()

# Concentration Series
with open(indir1 + 'pm_nat_conc_dict', 'rb') as f:
     pm_dict = cp.load(f)
f.close()

# Time Outside
with open(indir1 + 't_out_variation_net_benefits_no_race_1hr', 'rb') as f:
     time_nb = cp.load(f)
f.close()

# Adaptation Cost
with open(indir1 + 'ad_cost_variation_net_benefits_1hr', 'rb') as f:
     cost_nb = cp.load(f)
f.close()

# FINF
with open(indir1 + 'finf_variation_net_benefits_1hr', 'rb') as f:
     finf_nb = cp.load(f)
f.close()

#%%
# Define general info
clist = ['black','olivedrab', 'mediumorchid', 'grey']
fsize = 12
lsize = 12
xlabels_cost = ['-R,W\n$10.15\n$10.83','W\n      varies','R\n$32.62','+W\n$43.04','+R\n$75.65']
xticks_cost =[10.15, 17.76, 32.62, 43.04,75.65]
x_cost =[10.15,10.83, 17.76, 32.62, 43.04,75.65]
width = 0.2
space = 0.03
msize2 = 48
xticks_clim = [1,2]



# Iterate for each year:
for year in [2100]: #2050,
    # Set title
    title = 'fig_2_final_large'

    # Create climate variability plot and fill in data
    fig, ax = plt.subplots(2,2, figsize = (9, 6))
    
    # Create plot for climate variability
    ax[0,0].plot(pm_dict['REF'][year], climate_nb['REF'][year], color= clist[0], linestyle = '', marker='o', markersize=3, label = 'REF')
    ax[0,0].plot(pm_dict['P37'][year], climate_nb['P37'][year], color= clist[1], linestyle = '', marker='x', markersize=3, label = 'P3.7')
    
    ax[0,0].set_xticks([10, 11, 12, 13, 14, 15, 16])
    
    ax[0,0].set_xlabel(r'National Mean $PM_{2.5}$ Concentration $(\frac{\mu g}{m^3})$', fontsize = fsize)
    ax[0,0].set_title('Climate Variability', fontsize = fsize)
    ax[0,0].legend()
    
    ax[0,0].hlines(y = np.mean(climate_nb['REF'][year]), xmin=0,xmax =2200,linewidth=1, color=clist[0], linestyle = '--')
    ax[0,0].hlines(y = np.mean(climate_nb['P37'][year]), xmin=0,xmax =2200,linewidth=1, color=clist[1], linestyle = '--')
    
    ax[0,0].vlines(x = np.mean(pm_dict['REF'][year]), ymin=0,ymax =5000,linewidth=1, color=clist[0], linestyle = '--')
    ax[0,0].vlines(x = np.mean(pm_dict['P37'][year]), ymin=0,ymax =5000,linewidth=1, color=clist[1], linestyle = '--')
    
    
    if year == 2050:
        ax[0,0].set_ylim(0,750)
        ax[0,0].set_xlim(10,15)
        ax[0,0].annotate('(a)', xy = (12.5, 700), fontsize = fsize, horizontalalignment='center')
    if year == 2100:
        ax[0,0].set_ylim(0,1500)
        ax[0,0].set_xlim(10,16)
        ax[0,0].annotate('(a)', xy = (13.3, 1350), fontsize = fsize, horizontalalignment='center')
        
    # Create plot for time outside variability
    x1 = time_nb['REF'][year].keys()
    y1 = time_nb['REF'][year].values()
    x2 = time_nb['P37'][year].keys()
    y2 = time_nb['P37'][year].values()
    ax[0,1].plot(x1, y1, color= clist[0], marker='o', linestyle='', markersize=10, label = 'REF')
    ax[0,1].plot(x2, y2, color= clist[1], marker='x', linestyle='', markersize=10, label = 'P3.7')
    ax[0,1].set_xlabel('Hours', fontsize = fsize)
    ax[0,1].legend()
    ax[0,1].set_title('Outdoor Time Variation', fontsize = fsize)
    ax[0,1].hlines(y = time_nb['REF'][year][1], xmin = 0, xmax = 8, linewidth=1, color='grey', linestyle = '--')
    
    ax[0,1].set_xlim(0,8)
    ax[0,1].set_xticks([0, 2, 4, 6, 8])
    ax[0,1].set_xticklabels(['0', '2', '4', '6', '8'])
    
    if year == 2050:
        ax[0,1].set_ylim(0,4000)
        ax[0,1].vlines(x = 1,ymin = 0, ymax = 4000, linewidth=1, color='grey', linestyle = '--')
        ax[0,1].annotate('(b)', xy = (4, 3600), fontsize = fsize)
    if year == 2100:
        ax[0,1].set_ylim(0,7000)
        ax[0,1].vlines(x = 1,ymin = 0, ymax = 7000, linewidth=1, color='grey', linestyle = '--')
        ax[0,1].annotate('(b)', xy = (4, 6300), fontsize = fsize)
    
    # Create plot for adaptation cost variation
    x1 = cost_nb['REF'][year].keys()
    y1 = cost_nb['REF'][year].values()
    x2 = cost_nb['P37'][year].keys()
    y2 = cost_nb['P37'][year].values()
    ax[1,0].plot(x_cost, y1, color= clist[0], marker='o', linestyle='', markersize=10, label = 'REF')
    ax[1,0].plot(x_cost, y2, color= clist[1], marker='x', linestyle='', markersize=10, label = 'P3.7')
    ax[1,0].set_xlabel('Cost (USD 2020)', fontsize = fsize)
    ax[1,0].legend()
    ax[1,0].set_title('Cost Variation', fontsize = fsize)
    ax[1,0].set_xlim(0,80)
    ax[1,0].set_xticks(xticks_cost)
    ax[1,0].set_xticklabels(xlabels_cost)
    
    ax[1,0].hlines(y = cost_nb['REF'][year]['wage'], xmin = 0, xmax = 80, linewidth=1, color='grey', linestyle = '--')
    
    
    if year == 2050:
        ax[1,0].set_ylim(0,2000)
        ax[1,0].vlines(x = 17.76,ymin = 0, ymax = 3500, linewidth=1, color='grey', linestyle = '--')
        ax[1,0].annotate('(c)', xy = (38, 1800), fontsize = fsize)
    if year == 2100:
        ax[1,0].set_ylim(0,3000)
        ax[1,0].vlines(x = 17.76,ymin = 0, ymax = 5000, linewidth=1, color='grey', linestyle = '--')
        ax[1,0].annotate('(c)', xy = (38,2700), fontsize = fsize)
        
    # Create plot for FINF variation
    x1 = finf_nb['REF'][year].keys()
    y1 = finf_nb['REF'][year].values()
    x2 = finf_nb['P37'][year].keys()
    y2 = finf_nb['P37'][year].values()
    ax[1,1].plot(x1, y1, color= clist[0], marker='o', linestyle='', markersize=10, label = 'REF')
    ax[1,1].plot(x2, y2, color= clist[1], marker='x', linestyle='', markersize=10, label = 'P3.7')
    ax[1,1].set_xlabel('FINF', fontsize = fsize)
    ax[1,1].legend()
    ax[1,1].set_title('Infiltration Variation', fontsize = fsize)
    ax[1,1].hlines(y = finf_nb['REF'][year][0.2], xmin = 0, xmax = 0.8, linewidth=1, color='grey', linestyle = '--')
    ax[1,1].set_xlim(0,0.8)
    if year == 2050:
        ax[1,1].set_ylim(0,1000)
        ax[1,1].vlines(x = 0.2,ymin = 0, ymax = 2000, linewidth=1, color='grey', linestyle = '--')
        ax[1,1].annotate('(d)', xy = (0.4, 900), fontsize = fsize)
    if year == 2100:
        ax[1,1].set_ylim(0,1500)
        ax[1,1].vlines(x = 0.2,ymin = 0, ymax = 3000, linewidth=1, color='grey', linestyle = '--')
        ax[1,1].annotate('(d)', xy = (0.4, 1350), fontsize = fsize)
    
    plt.tight_layout()
    fig.text(0, 0.5, 'National Mean Per Capita Net Benefits of Adaptation (USD 2020)', va='center', rotation='vertical', fontsize = fsize)
    
    if year == 2050:
        fig.suptitle('2050', fontsize = 16)
    if year == 2100:
        fig.suptitle('2100', fontsize = 16)
        
    # plt.savefig(outdir + title + '.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')

