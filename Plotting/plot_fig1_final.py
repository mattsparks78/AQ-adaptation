# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:14:23 2022

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
indir1 = r".\Data\Rational Actor Infiltration\\"
indir2 = r".\Data\Future AQ\\"
indir3 = r".\Data\General\\"
indir4 = r".\Data\Social Learning Model\\"
indir5 = r".\Data\Historical Model PM\\"
indir6 = r'.\Data\cb_2018_us_nation_5m\\'
indir7 = r'.\Data\BenMap Grid\\'
outdir = r".\Figures\\"

# Set infiltration rate based on empirical data/for sensitvity analysis
inf_rate = 0.2

# Load general information
# Policy cases
pol_list = ['REF', 'P45', 'P37']

# IC cases
IC_list = ['IC' + str(i) for i in range(1, 6)]

# Population
pop = pd.read_csv(indir3 + 'pop_lepeule.csv')

# List of cells
rc = list(zip(pop['ROW'], pop['COL']))

# List of years
years = [range(2036, 2066), range(2086, 2116)]

# Population proportion in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (
        pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])

# Create dict to hold national pop-weighted aqi days by pol/pol_year
nat_aqi_days = {}

# Load AQI days pickles
for pol in pol_list:
    nat_aqi_days[pol] = {}
    filename = pol + '_tot_bad_2050_2100_d_pm_only'
    with open(indir2 + filename, 'rb') as f:
        aqi_all = cp.load(f)
    f.close()
    for pol_year in [2050, 2100]:
        nat_aqi_days[pol][pol_year] = []
        if pol_year == 2050:
            year_list = years[0]
        else:
            year_list = years[1]
        for year in year_list:
            year_tot = []
            for IC in IC_list:
                IC_tot = []
                for r, c in rc:
                    IC_tot.append(pop_prop[(r, c)] *
                                  aqi_all[str(year)][IC][(r, c)])
                year_tot.append(sum(IC_tot))
            nat_aqi_days[pol][pol_year].append(
                round(sum(year_tot) / len(year_tot)))


# Load AQI days and adaptation days for REF 2000
with open(indir5 + 'pm_conc_base_model_2000', 'rb') as f:
    pm_2000 = cp.load(f)
f.close()

AQ_high_2000 = {}
for r, c in rc:
    AQ_high_2000[(r, c)] = []
    for yr in range(1981, 2011):
        AQ_high_2000[(r, c)].append(
            len([x for x in pm_2000[yr][(r, c)] if x >= 35.5]))
    AQ_high_2000[(r, c)] = round(
        sum(AQ_high_2000[(r, c)]) / len(AQ_high_2000[(r, c)]))

pm_high_2000 = []

for yr in range(1981, 2011):
    bad_year = []
    for r, c in rc:
        bad = len([x for x in pm_2000[yr][(r, c)] if x >= 35.5])
        bad_year.append(bad * pop_prop[(r, c)])
    pm_high_2000.append(round(sum(bad_year)))

nat_AQI_2000 = sum(pm_high_2000) / len(pm_high_2000)
# %% Calculate contours for aqi days

with open(indir2 + 'grid_days_bad_2050_pm_100', 'rb') as f:
    aqi_gdf_2050 = cp.load(f)
f.close()

with open(indir2 + 'grid_days_bad_2100_pm_100', 'rb') as f:
    aqi_gdf_2100 = cp.load(f)
f.close()

# Convert REF 2000 to a dataframe
aqi_gdf_2000 = copy.deepcopy(aqi_gdf_2050[['ROW', 'COL']])
aqi_gdf_2000['REF Days'] = 0
for r, c in rc:
    aqi_gdf_2000.loc[(aqi_gdf_2000['ROW'] == r) & (
        aqi_gdf_2000['COL'] == c), 'REF Days'] = AQ_high_2000[(r, c)]

aqi_gdf_2050['diff'] = aqi_gdf_2050['REF Days'] - aqi_gdf_2000['REF Days']
aqi_gdf_2100['diff'] = aqi_gdf_2100['REF Days'] - aqi_gdf_2000['REF Days']

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
        crs='epsg:4326')


mymax = 50

norm = mpl.colors.Normalize(vmin=0, vmax=mymax)

us_bound = gpd.read_file(indir6 + 'cb_2018_us_nation_5m.shp')
us_bound = us_bound.to_crs(epsg=4326)

aqi_gdf_dict = {}
for year in [2050, 2100]:

    # Create array of x and y values
    lon_arr = aqi_gdf_2050['Lon'].unique()
    lat_arr = aqi_gdf_2050['Lat'].unique()

    # Change datatype of lon/lat arrays
    lon_arr = lon_arr.astype('float')
    lat_arr = lat_arr.astype('float')

    # Sort arrays
    lon_arr = np.sort(lon_arr)
    lat_arr = np.sort(lat_arr)

    # Create empty array to store pol ben values
    ben_arr = np.zeros((len(lat_arr), len(lon_arr)))

    if year == 2050:
        for i in range(len(lat_arr)):
            for j in range(len(lon_arr)):
                b = aqi_gdf_2050.loc[(aqi_gdf_2050['Lon'] == str(lon_arr[j])) & (
                    aqi_gdf_2050['Lat'] == str(lat_arr[i])), 'diff']
                if len(b) == 0:
                    pass
                else:
                    ben_arr[i, j] = aqi_gdf_2050.loc[(aqi_gdf_2050['Lon'] == str(lon_arr[j])) & (
                        aqi_gdf_2050['Lat'] == str(lat_arr[i])), 'diff'].item()

    if year == 2100:
        for i in range(len(lat_arr)):
            for j in range(len(lon_arr)):
                b = aqi_gdf_2100.loc[(aqi_gdf_2100['Lon'] == str(lon_arr[j])) & (
                    aqi_gdf_2100['Lat'] == str(lat_arr[i])), 'diff']
                if len(b) == 0:
                    pass
                else:
                    ben_arr[i, j] = aqi_gdf_2100.loc[(aqi_gdf_2100['Lon'] == str(lon_arr[j])) & (
                        aqi_gdf_2100['Lat'] == str(lat_arr[i])), 'diff'].item()

    lon_arr = sp.ndimage.zoom(lon_arr, 5)
    lat_arr = sp.ndimage.zoom(lat_arr, 5)
    ben_arr = sp.ndimage.zoom(ben_arr, 5)
    ben_arr[ben_arr < 0] = 0

    collec_poly = plt.contourf(
        lon_arr, lat_arr, ben_arr, levels=25, cmap='Reds', norm=norm)

    test = collec_to_gdf(collec_poly)

    ticks = [0, mymax*1/5, mymax*2/5, mymax*3/5, mymax*4/5, mymax]
    label_aqi = 'Extra Air Quality Alert Days Per Year\n Compared to Start-of-Century'

    colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
    test['RGBA'] = colors
    test = test.to_crs(epsg=4326)
    for i in range(len(test.geometry)):
        test['geometry'][i] = test['geometry'][i].intersection(
            us_bound['geometry'][0])

    aqi_gdf_dict[year] = test

# Plot regional boundaries
grid_reg = gpd.read_file(indir7 + 'US_census_regions.shp')

tick_labels = [str(0), str(round(mymax*1/5)), str(round(mymax*2/5)),
               str(round(mymax*3/5)), str(round(mymax*4/5)), str(round(mymax))]
# %%
# Calculate distribution of AQI days by race, hh income +- median, and proportion with/without AC
# Load grid cell racial breakdown
with open(indir3 + 'grid_race_dict', 'rb') as f:
    grid_race_dict = cp.load(f)
f.close()

# Load grid cell per capita income breakdown
with open(indir3 + 'grid_pc_inc_dict', 'rb') as f:
    grid_inc_dict = cp.load(f)
f.close()

# Combine others into all others
for r, c in rc:
    grid_race_dict[(r, c)]['All Others'] = grid_race_dict[(r, c)]['Asian'] + grid_race_dict[(r, c)]['Native American'] + \
        grid_race_dict[(r, c)]['Pacific'] + grid_race_dict[(r, c)
                                                           ]['Other'] + grid_race_dict[(r, c)]['Two Plus']

# Create income for all others
for r, c in rc:
    grid_inc_dict[(r, c)]['Income']['All Others'] = []
    for race in ['Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        race_pop_prop = grid_race_dict[(
            r, c)][race] / grid_race_dict[(r, c)]['All Others']
        race_inc = grid_inc_dict[(r, c)]['Income'][race]
        grid_inc_dict[(r, c)]['Income']['All Others'].append(
            race_inc * race_pop_prop)
    grid_inc_dict[(r, c)]['Income']['All Others'] = sum(
        grid_inc_dict[(r, c)]['Income']['All Others'])

# Create dict with grid/race breakdown
tot_pop_list = []
for r, c in rc:
    tot_pop_list.append(grid_race_dict[(r, c)]['Total'])
tot_pop = sum(tot_pop_list)

tot_race_pop = {}
for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    tot_race_pop[race] = []
    for r, c in rc:
        tot_race_pop[race].append(grid_race_dict[(r, c)][race])
    tot_race_pop[race] = sum(tot_race_pop[race])

temp_list = []
for race in ['Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    temp_list.append(tot_race_pop[race])
tot_race_pop['All Others'] = sum(temp_list)

grid_race_prop = {}
for r, c in rc:
    grid_race_prop[(r, c)] = {}
    for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
        grid_race_prop[(r, c)][race] = grid_race_dict[(r, c)
                                                      ][race] / tot_race_pop[race]

# Create df of income data
grid_inc_df = pd.DataFrame()
grid_inc_df['ROW'] = aqi_gdf_2050['ROW']
grid_inc_df['COL'] = aqi_gdf_2050['COL']

for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    grid_inc_df[race] = 0

for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    for r, c in rc:
        grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
            grid_inc_df['COL'] == c), race] = grid_inc_dict[(r, c)]['Income'][race]

# Calcuate average national median income from 2016 - 2020 (same timeframe as census data)
pc_inc_list = []

for race in ['White', 'Black', 'Asian', 'Native American', 'Pacific', 'Other', 'Two Plus']:
    for r, c in rc:
        pcinc = [grid_inc_dict[(r, c)]['Income'][race]]
        popu = grid_race_dict[(r, c)][race]
        pop_frac = max(round((popu / 1000)), 1)
        pc_inc_list.extend(pcinc*pop_frac)

pc_inc_list = np.array(pc_inc_list)

med_pc_inc = np.median(pc_inc_list)


# Compare cell hh income to national median
for race in ['White', 'Black', 'All Others']:
    grid_inc_df[race] = 0

for race in ['White', 'Black', 'All Others']:
    for r, c in rc:
        grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
            grid_inc_df['COL'] == c), race] = grid_inc_dict[(r, c)]['Income'][race]

for race in ['White', 'Black', 'All Others']:
    grid_inc_df['Comp' + race] = 0

for race in ['White', 'Black', 'All Others']:
    for r, c in rc:
        if grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() >= med_pc_inc:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'Above'
        elif 0 < grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() < med_pc_inc:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'Below'
        elif grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() == 0:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'None'


# Bring in PM concentrations for REF 2050 and 2100
with open(indir2 + 'REF_PM25_2050_2100_d', 'rb') as f:
    future_pm = cp.load(f)
f.close()

# Break down future pm concentrations by grid cell
future_pm_cell = {}
for year in [2050, 2100]:
    if year == 2050:
        year_list = range(2036, 2066)
    elif year == 2100:
        year_list = range(2086, 2116)
    future_pm_cell[year] = {}
    for r, c in rc:
        future_pm_cell[year][(r, c)] = []
        for yr in year_list:
            for IC in IC_list:
                future_pm_cell[year][(r, c)].append(
                    np.mean(future_pm[str(yr)][IC][(r, c)]))
        future_pm_cell[year][(r, c)] = sum(
            future_pm_cell[year][(r, c)]) / len(future_pm_cell[year][(r, c)])

# Calculate 150 year average national pm concentration
nat_pm = {}
for year in [2050, 2100]:
    nat_pm[year] = []
    for r, c in rc:
        nat_pm[year].append(future_pm_cell[year][(r, c)] * pop_prop[(r, c)])
    nat_pm[year] = sum(nat_pm[year])


# %%
# Prepare data for CDF Plot for AQI days by race
# Calc CDF for race distribution
cdf_conc_2050 = [0, 5, 10, 15, 20, 25, 30, 35]
cdf_conc_2100 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

race_cdf = {}
for year in [2050, 2100]:
    race_cdf[year] = {}
    if year == 2050:
        cdf_conc = cdf_conc_2050
    elif year == 2100:
        cdf_conc = cdf_conc_2100
    for race in ['White', 'Black', 'All Others']:
        race_cdf[year][race] = []
        for con in cdf_conc:
            p_race = 0
            for r, c in rc:
                if race in ['White', 'Black']:
                    popu = grid_race_dict[(r, c)][race]
                elif race == 'All Others':
                    popu = grid_race_dict[(r, c)]['Asian'] + grid_race_dict[(r, c)]['Native American']  \
                        + grid_race_dict[(r, c)]['Other'] + grid_race_dict[(r, c)]['Pacific']  \
                        + grid_race_dict[(r, c)]['Two Plus']
                if year == 2050:
                    days_future = aqi_gdf_2050.loc[(aqi_gdf_2050['ROW'] == r) & (
                        aqi_gdf_2050['COL'] == c), 'REF Days'].item()
                elif year == 2100:
                    days_future = aqi_gdf_2100.loc[(aqi_gdf_2100['ROW'] == r) & (
                        aqi_gdf_2100['COL'] == c), 'REF Days'].item()
                days_past = aqi_gdf_2000.loc[(aqi_gdf_2000['ROW'] == r) & (
                    aqi_gdf_2000['COL'] == c), 'REF Days'].item()
                if days_future - days_past >= con:
                    p_race += popu

            race_cdf[year][race].append(
                round((p_race / tot_race_pop[race]) * 100))

# %%
# Prepare data for CDF Plot of AQI days by +- national residential FINF
# Load pickle
with open(indir3 + 'grid_ach50_all', 'rb') as f:
    grid_ACH = cp.load(f)
f.close()

# Calculate national average
nat_ACH = []
for r, c in rc:
    nat_ACH.append(grid_ACH[(r, c)] * pop_prop[(r, c)])
nat_ACH = sum(nat_ACH)

# Create df of ACH data
grid_ach_df = pd.DataFrame()
grid_ach_df['ROW'] = grid_inc_df['ROW']
grid_ach_df['COL'] = grid_inc_df['COL']
grid_ach_df['Comp'] = 0

# Compare cell ACH to national mean
for r, c in rc:
    if grid_ACH[(r, c)] >= nat_ACH:
        grid_ach_df.loc[(grid_ach_df['ROW'] == r) & (
            grid_ach_df['COL'] == c), 'Comp'] = 'Above'
    elif grid_ACH[(r, c)] < nat_ACH:
        grid_ach_df.loc[(grid_ach_df['ROW'] == r) & (
            grid_ach_df['COL'] == c), 'Comp'] = 'Below'

# Plot CDF for AC/No AC
cdf_conc_2050 = [0, 5, 10, 15, 20, 25, 30, 35]
cdf_conc_2100 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

ACH_pop = {}
for comp in ['Above', 'Below']:
    ACH_pop[comp] = []
    s = grid_ach_df[grid_ach_df['Comp'] == comp]
    temp_rc = list(zip(s['ROW'], s['COL']))
    for r, c in temp_rc:
        ACH_pop[comp].append(grid_race_dict[(r, c)]['Total'])
    ACH_pop[comp] = sum(ACH_pop[comp])

ACH_cdf = {}
for comp in ['Above', 'Below']:
    ACH_cdf[comp] = {}
    s = grid_ach_df[grid_ach_df['Comp'] == comp]
    temp_rc = list(zip(s['ROW'], s['COL']))
    for year in [2050, 2100]:
        ACH_cdf[comp][year] = []
        if year == 2050:
            cdf_conc = cdf_conc_2050
        elif year == 2100:
            cdf_conc = cdf_conc_2100
        for con in cdf_conc:
            p_AC = 0
            for r, c in temp_rc:
                popu = grid_race_dict[(r, c)]['Total']
                if year == 2050:
                    days_future = aqi_gdf_2050.loc[(aqi_gdf_2050['ROW'] == r) & (
                        aqi_gdf_2050['COL'] == c), 'REF Days'].item()
                elif year == 2100:
                    days_future = aqi_gdf_2100.loc[(aqi_gdf_2100['ROW'] == r) & (
                        aqi_gdf_2100['COL'] == c), 'REF Days'].item()
                days_past = aqi_gdf_2000.loc[(aqi_gdf_2000['ROW'] == r) & (
                    aqi_gdf_2000['COL'] == c), 'REF Days'].item()
                if days_future - days_past >= con:
                    p_AC += popu

            ACH_cdf[comp][year].append(round((p_AC / ACH_pop[comp]) * 100))


# %%
# Prepare data for histogram plot of AQI days by += national median income
# Create column for overall income
grid_inc_df['All'] = 0
for r, c in rc:
    all_inc = []
    for race in ['White', 'Black', 'All Others']:
        prop_race = grid_race_dict[(r, c)][race] / \
            grid_race_dict[(r, c)]['Total']
        race_inc = race_inc = grid_inc_dict[(r, c)]['Income'][race]
        all_inc.append(prop_race * race_inc)
    grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
        grid_inc_df['COL'] == c), 'All'] = sum(all_inc)


# Compare cell hh income to national median
for race in ['White', 'Black', 'All Others', 'All']:
    grid_inc_df['Comp' + race] = 0

for race in ['White', 'Black', 'All Others', 'All']:
    for r, c in rc:
        if grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() >= med_pc_inc:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'Above'
        elif 0 < grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() < med_pc_inc:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'Below'
        elif grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (grid_inc_df['COL'] == c), race].item() == 0:
            grid_inc_df.loc[(grid_inc_df['ROW'] == r) & (
                grid_inc_df['COL'] == c), 'Comp' + race] = 'None'
# Create lists of AQI days for above/below median income
# Create lists of increase in AQI days for above/below
inc_days = {}
for comp in ['Above', 'Below']:
    inc_days[comp] = {}
    for year in [2050, 2100]:
        inc_days[comp][year] = []
        s = grid_inc_df[grid_inc_df['CompAll'] == comp]
        rc_sub = list(zip(s['ROW'], s['COL']))
        for r, c in rc_sub:
            if year == 2050:
                aqdiff = [aqi_gdf_2050.loc[(aqi_gdf_2050['ROW'] == r) & (aqi_gdf_2050['COL'] == c), 'REF Days'].item() -
                          aqi_gdf_2000.loc[(aqi_gdf_2000['ROW'] == r) & (aqi_gdf_2000['COL'] == c), 'REF Days'].item()]
            elif year == 2100:
                aqdiff = [aqi_gdf_2100.loc[(aqi_gdf_2100['ROW'] == r) & (aqi_gdf_2100['COL'] == c), 'REF Days'].item() -
                          aqi_gdf_2000.loc[(aqi_gdf_2000['ROW'] == r) & (aqi_gdf_2000['COL'] == c), 'REF Days'].item()]
            popu = grid_race_dict[(r, c)]['Total']
            pop_frac = max(round((popu / 1000)), 1)
            inc_days[comp][year].extend(aqdiff*pop_frac)

# Convert to np arrays
for comp in ['Above', 'Below']:
    for year in [2050, 2100]:
        inc_days[comp][year] = np.array(inc_days[comp][year])

# %%
# Create figure with three subfigures
fig = plt.figure(figsize=(12, 6))

subfigs = fig.subfigures(1, 3, wspace=0, width_ratios=[1, 1, 1])

axsLeft = subfigs[0].subplots(2, 1)
axsMid = subfigs[1].subplots(2, 1)
axsRight = subfigs[2].subplots(2, 1)
plt.subplots_adjust(hspace=0.25, wspace=0.25)
# fig.tight_layout()

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)


# Plot with broken y axis
labels = ['2050', '2100']
width = 0.25
space = 0.2
xticks = [2, 4, 6]
clist = ['black', (228/255, 180/255, 41/255), (127/255, 127/255, 127/255)]
xticklabels = ['2000', '2050', '2100']
lsize = 12

bgcolors = ['#FFFFFF', (201/255, 231/255, 255/255),
            (201/255, 231/255, 255/255)]

# Create first plot
axsLeft[0].boxplot(pm_high_2000, positions=[xticks[0]], widths=width, showfliers=False,
                   patch_artist=True, boxprops=dict(facecolor=clist[0], color=clist[0]), capprops=dict(color=clist[0]),
                   whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
# REF
r = axsLeft[0].boxplot(nat_aqi_days['REF'][2050], positions=[xticks[1]], widths=width, showfliers=False,
                       patch_artist=True, boxprops=dict(facecolor=clist[0], color=clist[0]), capprops=dict(color=clist[0]),
                       whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
axsLeft[0].boxplot(nat_aqi_days['REF'][2100], positions=[xticks[2]], widths=width, showfliers=False,
                   patch_artist=True, boxprops=dict(facecolor=clist[0], color=clist[0]), capprops=dict(color=clist[0]),
                   whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))

axsLeft[0].set_ylabel('Air Quality Alert Days per Year', fontsize=lsize)
# axsLeft[0].set_title('Air Quality\n Advisories', fontsize = lsize)
axsLeft[0].set_ylim([0, 30])
axsLeft[0].set_yticks([0, 10, 20, 30])
axsLeft[0].set_yticklabels(['0', '10', '20', '30'])
# axsLeft[0].legend([r["boxes"][0]], ['REF'], loc='upper left', fontsize=12)
axsLeft[0].set_xlim([0, 8])
axsLeft[0].set_xticks(xticks)
axsLeft[0].set_xticklabels(xticklabels)
axsLeft[0].set_xlabel('Year', fontsize=lsize)

axsLeft[0].annotate('(a)', (3.7, 26), fontsize=lsize)

# Plot bottom of left subfigure. Historam of AQI days by income
axsLeft[1].hist(inc_days['Above'][2100], bins=39,
                stacked=True, density=True, color='darkgreen')
axsLeft[1].hist(inc_days['Below'][2100], bins=45,
                stacked=True, density=True, color='palegreen')
axsLeft[1].set_xlim(0, 45)

axsLeft[1].legend({'Above': "darkgreen", 'Below': "palegreen"}, fontsize=lsize)

axsLeft[1].set_ylabel('Portion of population', fontsize=lsize)
axsLeft[1].annotate('Average Above Median Income = ' + str(round(np.mean(inc_days['Above'][2100]))) + '\n'
                    'Average Below Median Income = ' + str(round(np.mean(inc_days['Below'][2100]))), (5, 0.2))
axsLeft[1].set_xlabel('Extra Air Quality Alert Days Per Year\n Compared to Start-of-Century', fontsize=lsize)
axsLeft[1].annotate('(b)', (20, 0.35), fontsize=lsize)

# Plot second subfigure
aqi_gdf_dict[2050].plot(color=aqi_gdf_dict[2050]['RGBA'], ax=axsMid[0])
aqi_gdf_dict[2100].plot(color=aqi_gdf_dict[2100]['RGBA'], ax=axsMid[1])

axsMid[0].set_title('2050', fontsize=lsize)
axsMid[1].set_title('2100', fontsize=lsize)

axsMid[0].annotate('(c)', (-126, 23.5), fontsize=lsize)
axsMid[1].annotate('(d)', (-126, 23.5), fontsize=lsize)

axsMid[0].set_axis_off()
axsMid[1].set_axis_off()

grid_reg.boundary.plot(ax=axsMid[0], color='black')
grid_reg.boundary.plot(ax=axsMid[1], color='black')

ax4 = subfigs[1].add_axes([0.05, .075, .9, 0.03])
cbar = mpl.colorbar.ColorbarBase(
    ax4, cmap='Reds', orientation='horizontal', ticks=ticks, norm=norm)
cbar.set_label(label_aqi, fontsize=lsize)
cbar.ax.set_xticklabels(tick_labels)


# Plot third subfigure
# CDF of s by race
# 2050
a, = axsRight[0].plot(cdf_conc_2050, race_cdf[2050]['White'],
                      color='blue', marker='o', label='White')
a2, = axsRight[0].plot(cdf_conc_2050, race_cdf[2050]
                       ['White'], color='blue', marker='o', label='2050')
b, = axsRight[0].plot(cdf_conc_2050, race_cdf[2050]['Black'],
                      color='red', marker='x', label='Black')
c, = axsRight[0].plot(cdf_conc_2050, race_cdf[2050]['All Others'],
                      color='lawngreen', marker='^', label='All Others')

axsRight[0].set_ylim(0, 100)
axsRight[0].set_xlim(0, 45)
axsRight[0].set_ylabel('Percent of population', fontsize=lsize)
# 2100
d, = axsRight[0].plot(cdf_conc_2100, race_cdf[2100]['White'],
                      color='darkblue', marker='o', linestyle='--', label='2100')
axsRight[0].plot(cdf_conc_2100, race_cdf[2100]['Black'],
                 color='darkred', marker='x', linestyle='--')
axsRight[0].plot(cdf_conc_2100, race_cdf[2100]['All Others'],
                 color='darkolivegreen', marker='^', linestyle='--')


first_legend = axsRight[0].legend(handles=[a, b, c], fontsize=lsize)
axsRight[0].add_artist(first_legend)

axsRight[0].legend(handles=[a2, d], loc=(0.16, 0.67), fontsize=lsize)

# CDF of s by FINF
# 2050
axsRight[1].plot(cdf_conc_2050, ACH_cdf['Above']
                 [2050], color='darkgrey', marker='*')
axsRight[1].plot(cdf_conc_2050, ACH_cdf['Below']
                 [2050], color='cyan', marker='s')

axsRight[1].set_ylim(0, 100)
axsRight[1].set_xlim(0, 45)


# 2100
axsRight[1].plot(cdf_conc_2100, ACH_cdf['Above'][2100],
                 color='dimgrey', marker='*', linestyle='--')
axsRight[1].plot(cdf_conc_2100, ACH_cdf['Below'][2100],
                 color='darkcyan', marker='s', linestyle='--')
axsRight[1].legend({'Leakier 2050': "darkgrey", 'Tighter 2050': "cyan",
                   'Leakier 2100': "dimgrey", 'Tighter 2100': "darkcyan"}, fontsize=lsize)
axsRight[1].set_xlabel(
    'Extra Air Quality Alert Days Per Year\n Compared to Start-of-Century', fontsize=lsize)
axsRight[1].set_ylabel('Percent of population', fontsize=lsize)

axsRight[0].annotate('(e)', (3, 10), fontsize=lsize)
axsRight[1].annotate('(f)', (3, 10), fontsize=lsize)

title = 'fig_1_final'

# plt.savefig(outdir + title + '.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')
