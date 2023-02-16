# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 09:01:45 2022

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
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"

# Load Future PM data
with open(indir1 + 'REF_PM25_2050_2100_d', 'rb') as f:
     pm_REF = cp.load(f)
f.close()

with open(indir1 + 'P45_PM25_2050_2100_d', 'rb') as f:
     pm_P45 = cp.load(f)
f.close()

with open(indir1 + 'P37_PM25_2050_2100_d', 'rb') as f:
     pm_P37 = cp.load(f)
f.close()

# Calculate proportion of AQI days over 100 per year
aqi_REF = {}
aqi_P45 = {}
aqi_P37 = {}

for Yr in pm_P37.keys():
    aqi_REF[Yr] = {}
    aqi_P45[Yr] = {}
    aqi_P37[Yr] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        aqi_REF[Yr][IC] = {}
        aqi_P45[Yr][IC] = {}
        aqi_P37[Yr][IC] = {}
        for cell in pm_REF[Yr][IC].keys():
            aqi_REF[Yr][IC][cell] = []
            aqi_P45[Yr][IC][cell] = []
            aqi_P37[Yr][IC][cell] = []
            for i in range(len(pm_REF[Yr][IC][cell])):
                con_pm_REF = pm_REF[Yr][IC][cell][i]
                con_pm_P45 = pm_P45[Yr][IC][cell][i]
                con_pm_P37 = pm_P37[Yr][IC][cell][i]
                
                #REF
                if con_pm_REF <= 12:
                    aqi_pm_REF = round((4.17*con_pm_REF))
                elif 12 < con_pm_REF <= 35.4:
                    aqi_pm_REF = round((2.1*(con_pm_REF - 12) + 51))
                elif 35.4 < con_pm_REF <= 55.4:
                    aqi_pm_REF = round((2.46*(con_pm_REF - 35.4) + 101))
                elif 55.4 < con_pm_REF <= 150.4:
                    aqi_pm_REF = round((0.52*(con_pm_REF-55.4) + 151))
                elif 150.4 < con_pm_REF <= 250.4:
                    aqi_pm_REF = round((0.99*(con_pm_REF-150.4) + 201))
                else:
                    aqi_pm_REF = round((0.8*(con_pm_REF-250.4) + 301))

                #P45
                if con_pm_P45 <= 12:
                    aqi_pm_P45 = round((4.17*con_pm_P45))
                elif 12 < con_pm_P45 <= 35.4:
                    aqi_pm_P45 = round((2.1*(con_pm_P45 - 12) + 51))
                elif 35.4 < con_pm_P45 <= 55.4:
                    aqi_pm_P45 = round((2.46*(con_pm_P45 - 35.4) + 101))
                elif 55.4 < con_pm_P45 <= 150.4:
                    aqi_pm_P45 = round((0.52*(con_pm_P45-55.4) + 151))
                elif 150.4 < con_pm_P45 <= 250.4:
                    aqi_pm_P45 = round((0.99*(con_pm_P45-150.4) + 201))
                else:
                    aqi_pm_P45 = round((0.8*(con_pm_P45-250.4) + 301))

                #P37
                if con_pm_P37 <= 12:
                    aqi_pm_P37 = round((4.17*con_pm_P37))
                elif 12 < con_pm_P37 <= 35.4:
                    aqi_pm_P37 = round((2.1*(con_pm_P37 - 12) + 51))
                elif 35.4 < con_pm_P37 <= 55.4:
                    aqi_pm_P37 = round((2.46*(con_pm_P37 - 35.4) + 101))
                elif 55.4 < con_pm_P37 <= 150.4:
                    aqi_pm_P37 = round((0.52*(con_pm_P37-55.4) + 151))
                elif 150.4 < con_pm_P37 <= 250.4:
                    aqi_pm_P37 = round((0.99*(con_pm_P37-150.4) + 201))
                else:
                    aqi_pm_P37 = round((0.8*(con_pm_P37-250.4) + 301))
                
                aqi_REF[Yr][IC][cell].append(aqi_pm_REF)
                aqi_P45[Yr][IC][cell].append(aqi_pm_P45)
                aqi_P37[Yr][IC][cell].append(aqi_pm_P37)
            
# Calculate proportion of bad air days in each scenario/cell
aqi_prop_REF = {}
aqi_prop_P45 = {}
aqi_prop_P37 = {}

for Yr in pm_P37.keys():
    aqi_prop_REF[Yr] = {}
    aqi_prop_P45[Yr] = {}
    aqi_prop_P37[Yr] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        aqi_prop_REF[Yr][IC] = {}
        aqi_prop_P45[Yr][IC] = {}
        aqi_prop_P37[Yr][IC] = {}
        for cell in pm_REF[Yr][IC].keys():
            aqi_prop_REF[Yr][IC][cell] = len([i for i in aqi_REF[Yr][IC][cell] if i > 100 ]) / 365
            aqi_prop_P45[Yr][IC][cell] = len([i for i in aqi_P45[Yr][IC][cell] if i > 100 ]) / 365
            aqi_prop_P37[Yr][IC][cell] = len([i for i in aqi_P37[Yr][IC][cell] if i > 100 ]) / 365
 
# Calculate proportion of bad air days in each scenario/cell
aqi_tot_REF = {}
aqi_tot_P45 = {}
aqi_tot_P37 = {}

for Yr in pm_P37.keys():
    aqi_tot_REF[Yr] = {}
    aqi_tot_P45[Yr] = {}
    aqi_tot_P37[Yr] = {}
    for IC in ['IC' + str(i) for i in range(1,6)]:
        aqi_tot_REF[Yr][IC] = {}
        aqi_tot_P45[Yr][IC] = {}
        aqi_tot_P37[Yr][IC] = {}
        for cell in pm_REF[Yr][IC].keys():
            aqi_tot_REF[Yr][IC][cell] = len([i for i in aqi_REF[Yr][IC][cell] if i > 100 ])
            aqi_tot_P45[Yr][IC][cell] = len([i for i in aqi_P45[Yr][IC][cell] if i > 100 ])
            aqi_tot_P37[Yr][IC][cell] = len([i for i in aqi_P37[Yr][IC][cell] if i > 100 ])
            



               
# Output pickles
with open(outdir + 'REF_AQI_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_REF,f)
f.close()

with open(outdir + 'P45_AQI_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_P45,f)
f.close()

with open(outdir + 'P37_AQI_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_P37,f)
f.close()

with open(outdir + 'REF_prop_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_prop_REF,f)
f.close()

with open(outdir + 'P45_prop_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_prop_P45,f)
f.close()

with open(outdir + 'P37_prop_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_prop_P37,f)
f.close()

with open(outdir + 'REF_tot_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_tot_REF,f)
f.close()

with open(outdir + 'P45_tot_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_tot_P45,f)
f.close()

with open(outdir + 'P37_tot_bad_2050_2100_d_pm_only', 'wb') as f:
      cp.dump(aqi_tot_P37,f)
f.close()


