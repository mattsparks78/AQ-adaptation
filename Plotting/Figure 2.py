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


# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir4 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Social Learning Model\\"
indir5 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor IF 0.5\\"
indir6 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor IF 0.25\\"
indir7 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Historical Model PM\\"
indir8 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\cb_2018_us_nation_5m\\'
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\Rational\\"

# Load general information
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

# Population proportion in each cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])

#Create dict to hold national pop-weighted aqi days by pol/pol_year
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
                    IC_tot.append(pop_prop[(r, c)] * aqi_all[str(year)][IC][(r, c)])    
                year_tot.append(sum(IC_tot))
            nat_aqi_days[pol][pol_year].append(round(sum(year_tot) / len(year_tot)))
        

# Load adaptation days pickles
nat_ad_days = {}
for adap_effic in [1, 0.5, 0.25]:
    nat_ad_days[adap_effic] = {}
    if adap_effic == 1:
        indir = indir1
    elif adap_effic == 0.5:
        indir = indir5
    else:
        indir = indir6
    for pol in pol_list:
        nat_ad_days[adap_effic][pol] = {}
        filename = 'adapt_count_all_' + pol
        with open(indir + filename, 'rb') as f:
             ad_days = cp.load(f)
        f.close()
        for pol_year in [2050, 2100]:
            nat_ad_days[adap_effic][pol][pol_year] = []
            if pol_year == 2050:
                year_list = years[0]
            else:
                year_list = years[1]
            for year in year_list:
                year_tot = []
                for IC in IC_list:
                    IC_tot = []
                    for r, c in rc:
                        IC_tot.append(pop_prop[(r, c)] * ad_days[year][IC][(r, c)])    
                    year_tot.append(sum(IC_tot))
                nat_ad_days[adap_effic][pol][pol_year].append(round(sum(year_tot) / len(year_tot)))
        
# Load AQI days and adaptation days for REF 2000
with open(indir7 + 'pm_conc_base_model_2000', 'rb') as f:
     pm_2000 = cp.load(f)
f.close()

pm_high_2000 = []

for yr in range(1981, 2011):
    bad_year = []
    for r, c in rc:
        bad = len([x for x in pm_2000[yr][(r, c)] if x >= 35.5 ])
        bad_year.append(bad * pop_prop[(r, c)])
    pm_high_2000.append(round(sum(bad_year)))


# Load adaptation count
with open(indir1 + 'adapt_count_all_REF_2000', 'rb') as f:
     adapt_2000 = cp.load(f)
f.close()

ad_2000 = []

for yr in range(1981, 2011):
    ad_year = []
    for r, c in rc:
        ad = adapt_2000[(r, c)][yr]
        ad_year.append(ad * pop_prop[(r, c)])
    ad_2000.append(round(sum(ad_year)))   
#%% Calculate contours for aqi days

with open(indir2 + 'grid_days_bad_2050_pm_100', 'rb') as f:
     aqi_gdf_2050 = cp.load(f)
f.close()

with open(indir2 + 'grid_days_bad_2100_pm_100', 'rb') as f:
     aqi_gdf_2100 = cp.load(f)
f.close()

aqi_gdf_2050['diff'] = aqi_gdf_2050['REF Days'] - aqi_gdf_2050['P37 Days']
aqi_gdf_2100['diff'] = aqi_gdf_2100['REF Days'] - aqi_gdf_2100['P37 Days']

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
 
    
mymax = 25

norm = mpl.colors.Normalize(vmin= 0, vmax= mymax)

us_bound = gpd.read_file(indir8 + 'cb_2018_us_nation_5m.shp')
us_bound = us_bound.to_crs(epsg=4326)

aqi_gdf_dict = {}
for year in [2050, 2100]:
    
    # Create array of x and y values
    lon_arr = aqi_gdf_2050['Lon'].unique()
    lat_arr = aqi_gdf_2050['Lat'].unique()
    
    # Change datatype of lon/lat arrays
    lon_arr = lon_arr.astype('float')
    lat_arr = lat_arr.astype('float')
    
    #Sort arrays
    lon_arr = np.sort(lon_arr)
    lat_arr = np.sort(lat_arr)
    
    # Create empty array to store pol ben values
    ben_arr = np.zeros((len(lat_arr), len(lon_arr)))
    
    if year == 2050:
        for i in range(len(lat_arr)):
            for j in range(len(lon_arr)):
                b = aqi_gdf_2050.loc[(aqi_gdf_2050['Lon'] == str(lon_arr[j])) & (aqi_gdf_2050['Lat'] == str(lat_arr[i])), 'diff']
                if len(b) == 0:
                    pass 
                else:
                    ben_arr[i, j] = aqi_gdf_2050.loc[(aqi_gdf_2050['Lon'] == str(lon_arr[j])) & (aqi_gdf_2050['Lat'] == str(lat_arr[i])), 'diff'].item() 
    
    if year == 2100:
        for i in range(len(lat_arr)):
            for j in range(len(lon_arr)):
                b = aqi_gdf_2100.loc[(aqi_gdf_2100['Lon'] == str(lon_arr[j])) & (aqi_gdf_2100['Lat'] == str(lat_arr[i])), 'diff']
                if len(b) == 0:
                    pass 
                else:
                    ben_arr[i, j] = aqi_gdf_2100.loc[(aqi_gdf_2100['Lon'] == str(lon_arr[j])) & (aqi_gdf_2100['Lat'] == str(lat_arr[i])), 'diff'].item()
                    
                    
    lon_arr = sp.ndimage.zoom(lon_arr, 5)
    lat_arr = sp.ndimage.zoom(lat_arr, 5)
    ben_arr = sp.ndimage.zoom(ben_arr, 5)
    ben_arr[ben_arr < 0] = 0
    
            
    collec_poly = plt.contourf(lon_arr, lat_arr, ben_arr, levels = 25, cmap = 'Oranges', norm = norm);
    
    test = collec_to_gdf(collec_poly)
    
    ticks=[0, mymax*1/5,mymax*2/5,mymax*3/5,mymax*4/5,mymax]
    label_aqi = 'Annual difference of\n Air QualityAdvisory days\n in REF due to '+  r'$PM_{2.5}$'    
    
    colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
    test['RGBA'] = colors
    test = test.to_crs(epsg=4326)
    for i in range(len(test.geometry)):
        test['geometry'][i] = test['geometry'][i].intersection(us_bound['geometry'][0])
    
    aqi_gdf_dict[year] = test



#%%
# Calculate contours for adapt days


# Set directories
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
indir4 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\US Boundary\\"
indir5 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\cb_2018_us_nation_5m\\'
indir6 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\\"

# Load future concentration data. REF, P3.7 or P4.5
with open(indir2 + 'REF_PM25_2050_2100_d', 'rb') as f:
    pm_out_REF = cp.load(f)
f.close()

with open(indir2 + 'P45_PM25_2050_2100_d', 'rb') as f:
    pm_out_P45 = cp.load(f)
f.close()

with open(indir2 + 'P37_PM25_2050_2100_d', 'rb') as f:
    pm_out_P37 = cp.load(f)
f.close()


# Define lists of years for 2050 and 2100
year_set = [list(range(2036, 2066)), list(range(2086, 2116))]

# Define list of Initial Conditions
IC_list = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5']

# Import hourly wage data
wage = pd.read_csv(indir1 + "wage_per_grid_all_b.csv")
# Only keep row, col, and year 2000 data
wage = wage[['ROW', 'COL', '2000']]
# Drop grid cells that aren't included in BenMap analysis
wage.drop([1, 219, 220], axis=0, inplace=True)
# Reset index back to starting at 0
wage.reset_index(inplace=True)

# Import population per grid cell
pop = pd.read_csv(indir1 + "pop_lepeule.csv")

# Create list of all grid cells
rc = list(zip(list(wage['ROW']), list(wage['COL'])))

# Calculate proportion of population in each grid cell
pop_prop = {}
for r, c in rc:
    pop_prop[(r, c)] = pop.loc[(pop['ROW'] == r) & (pop['COL'] == c), '2000'].item() / np.sum(pop['2000'])
    

# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')


# Set adaptation time in hours. Assuming full adaptation because it maximizes utility.
AT = 2

# Set relative risk
RR = 1.14

# Calculate attributable fraction
AF = (RR - 1) / RR

# Set VSL
VSL = 7400000

# Proportion of population who adapt. Assumed 1 here. Can update with replicator numbers.
ppa = 1


# Calculate adaptation days for all scenarios
# Initialize dictionaries
pm_out_base = {}
pm_out_pol = {}
adapt_count = {}

us_bound = gpd.read_file(indir5 + 'cb_2018_us_nation_5m.shp')
us_bound = us_bound.to_crs(epsg=4326)

# Iterate over both sets of years (for 2050 and 2100)
for pol in ['P37','P45']:
    adapt_count[pol] = {}
    for years in year_set:
        # Iterate over all years in year set
        if np.mean(years) < 2070:
            # Economic projection factor
            EPF = 2.01
            # Population projection factor
            PPF = 1.48
            # PM Mortality Factor (for baseline mortaltiy)
            PMMF = 1.12
            # Set VSL
            VSL = 7400000
            # Update VSL based on EPF
            VSL = VSL * EPF
        else:
            EPF = 3.94
            PPF = 1.73
            PMMF = 0.92
            # Set VSL
            VSL = 7400000
            # Update VSL based on EPF
            VSL = VSL * EPF
        for Yr in years:
            # Initialize next level of dictionaries
            pm_out_base[Yr] = {}
            pm_out_pol[Yr] = {}
            adapt_count[pol][Yr] = {}
            
    
            # Iterate over all initial conditions in each year
            for IC in IC_list:
                # Initialize next level of dictionaries
                pm_out_base[Yr][IC] = {}
                pm_out_pol[Yr][IC] = {}
                adapt_count[pol][Yr][IC] = {}
    
                # Iterate over all cells in each scenario
                for r, c in rc:
                    # Calculate pol mean outdoor PM concentration
                    if pol == 'P37':
                        pm_target = np.mean(pm_out_P37[str(Yr)][IC][(r, c)])
                    elif pol == 'P45':
                        pm_target = np.mean(pm_out_P45[str(Yr)][IC][(r, c)])
                    
                    
                    # Calculate baseline mean outdoor PM concentration
                    pm_o = np.mean(pm_out_REF[str(Yr)][IC][(r, c)])
                    # Store baseline mean outdoor PM concentration in dictionary
                    pm_out_base[Yr][IC][(r, c)] = pm_o
                    # Sort PM concentrations (Python is only lowest to highest :( )
                    pm_sort = np.sort(pm_out_REF[str(Yr)][IC][(r, c)])
                    # So we reverse the order of the array elements
                    pm_descend = copy.deepcopy(pm_sort[::-1])
                    # Calculate dPM threshold for adaptation
                    dPM_thresh = pm_o - pm_target
                    # Set adaptation count
                    count = 0
                    # Create deepcopy of pm concentrations in descending order
                    pm_test = copy.deepcopy(pm_descend)
                    # Creat dPM_total
                    dPM_total = 0
    
                    # Calculate the potential reduction of mean outdoor PM by adapting each day
                    for i in range(len(pm_descend)):
                        pm_test[i] = 0
                        pm_test_mean = np.mean(pm_test)
                        dPM_test = np.mean(pm_descend) - pm_test_mean
                        dPM_total += dPM_test
                        
                        if dPM_total < dPM_thresh:
                            pm_descend[i] = 0
                            count += 1
                        else:
                            break
    
                    # Output adaptation days count
                    adapt_count[pol][Yr][IC][(r, c)] = count
    
#%%
# Convert adapt count by cell to national pop-weighted average
adapt_nat = {}


adapt_cell_2050 = {}
adapt_cell_2100 = {}

for pol in ['P37', 'P45']:
    adapt_cell_2050[pol] = {}
    adapt_cell_2100[pol] = {}
    for r, c in rc:
        adapt_cell_2050[pol][(r, c)] = []
        adapt_cell_2100[pol][(r, c)] = []
        for years in year_set:
            for Yr in years:
                for IC in IC_list:
                    if Yr < 2070:
                        adapt_cell_2050[pol][(r, c)].append(adapt_count[pol][Yr][IC][(r, c)])
                    else:
                        adapt_cell_2100[pol][(r, c)].append(adapt_count[pol][Yr][IC][(r, c)])
                        
# Calculate range of national days
for year in [2050, 2100]:
    adapt_nat[year] = {}
    for pol in ['P37', 'P45']:
        adapt_nat[year][pol] = []
        if year == 2050:
            for i in range(0, 150):
                year_ad = []
                for r, c in rc:
                    prop_ref = pop_prop[(r, c)]
                    ad_cell = adapt_cell_2050[pol][(r, c)][i]
                    year_ad.append(prop_ref * ad_cell)
                adapt_nat[year][pol].append(round(sum(year_ad)))    
        else:
            for i in range(0, 150):
                year_ad = []
                for r, c in rc:
                    prop_ref = pop_prop[(r, c)]
                    ad_cell = adapt_cell_2100[pol][(r, c)][i]
                    year_ad.append(prop_ref * ad_cell)
                adapt_nat[year][pol].append(round(sum(year_ad))) 
            
mean_adapt_nat = {}
for year in [2050,2100]:
    mean_adapt_nat[year] = {}
    for pol in ['P37','P45']:
        mean_adapt_nat[year][pol] = round(np.mean(adapt_nat[year][pol]))
       


mean_adapt_cell_2050 = {}
mean_adapt_cell_2100 = {}
for pol in ['P37', 'P45']:
    mean_adapt_cell_2050[pol] = {}
    mean_adapt_cell_2100[pol] = {}
    for r, c in rc:
        mean_adapt_cell_2050[pol][(r, c)] = round(np.mean(adapt_cell_2050[pol][(r,c)]))
        mean_adapt_cell_2100[pol][(r, c)] = round(np.mean(adapt_cell_2100[pol][(r,c)]))


        
#%% Create contour plots for RA and SL

# Plot Rational Actor Map
# Import grid  gdf
grid = gpd.read_file(indir6 + "BenMapGridPoints.csv")
grid = grid.drop('geometry', axis = 1)
grid['geom'] = grid['geom'].apply(wkt.loads)
grid = grid.set_geometry('geom')
grid = grid.set_crs(epsg=4326)
grid['ROW'] = grid['ROW'].astype('int')
grid['COL'] = grid['COL'].astype('int')
grid = grid.drop([1, 219, 220])

# Plot adaptation days by grid cell
days_gdf = copy.deepcopy(grid)
for pol in ['P37', 'P45']:
    for year in [2050, 2100]:
        key = pol + ' ' + str(year)
        days_gdf[key] = 0

for pol in ['P37', 'P45']:
    for year in [2050, 2100]:
        key = pol + ' ' + str(year)
        if year == 2050:
            for r, c in rc:
                days_gdf.loc[(days_gdf['ROW'] == r) & (days_gdf['COL'] == c), key] = mean_adapt_cell_2050[pol][(r, c)]
        else:
            for r, c in rc:
                days_gdf.loc[(days_gdf['ROW'] == r) & (days_gdf['COL'] == c), key] = mean_adapt_cell_2100[pol][(r, c)] 






mymax = 25

norm = mpl.colors.Normalize(vmin= 0, vmax= mymax)

# Plot Contours


gdf_dict = {}
for pol in ['P37', 'P45']:
    gdf_dict[pol] = {}
    for year in [2050, 2100]:
        gdf_dict[pol][year] = {}
        key = pol + ' ' + str(year)
        
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
                b = days_gdf.loc[(days_gdf['Lon'] == str(lon_arr[j])) & (days_gdf['Lat'] == str(lat_arr[i])), key]
                if len(b) == 0:
                    pass 
                else:
                    ben_arr[i, j] = days_gdf.loc[(days_gdf['Lon'] == str(lon_arr[j])) & (days_gdf['Lat'] == str(lat_arr[i])), key].item() 
        
        lon_arr = sp.ndimage.zoom(lon_arr, 5)
        lat_arr = sp.ndimage.zoom(lat_arr, 5)
        ben_arr = sp.ndimage.zoom(ben_arr, 5)
        
                
        collec_poly = plt.contourf(lon_arr, lat_arr, ben_arr, levels = 25, cmap = 'Purples', norm = norm);
        
        test = collec_to_gdf(collec_poly)
        
        ticks=[0, mymax*1/5,mymax*2/5,mymax*3/5,mymax*4/5,mymax]
        label = 'Annual number of \nadaptation days to reach POL\n ' + r'$PM_{2.5}$' ' exposures under REF'
        tick_labels = [str(0), str(round(mymax*1/5)),str(round(mymax*2/5)),str(round(mymax*3/5)),str(round(mymax*4/5)),str(round(mymax))]
        
        
        colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
        test['RGBA'] = colors
        test = test.to_crs(epsg=4326)
        for i in range(len(test.geometry)):
            test['geometry'][i] = test['geometry'][i].intersection(us_bound['geometry'][0])
        
        gdf_dict[pol][year] = test










#%%
# Create figure with three subfigures
fig = plt.figure(figsize=(12, 6))

subfigs = fig.subfigures(1, 3, wspace=0)

axsLeft = subfigs[0].subplots(2, 2, gridspec_kw = {'height_ratios' : [1, 0.05]})
axsMid = subfigs[1].subplots(2, 1)
axsRight = subfigs[2].subplots(2, 1)
plt.subplots_adjust(hspace=0, wspace = 0.25)
fig.tight_layout()


### Plot with broken y axis
labels = ['2050', '2100']
width = 0.25
xticks = [2, 4, 6]
clist = ['black', 'red', 'blue']
xticklabels = ['2000', '2050', '2100']
lsize = 16


### Create first plot
axsLeft[0, 0].boxplot(pm_high_2000, positions = [xticks[0]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
# REF
r = axsLeft[0, 0].boxplot(nat_aqi_days['REF'][2050], positions = [xticks[1] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
axsLeft[0, 0].boxplot(nat_aqi_days['REF'][2100], positions = [xticks[2] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))

# P45
p_45 = axsLeft[0, 0].boxplot(nat_aqi_days['P45'][2050], positions = [xticks[1]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
axsLeft[0, 0].boxplot(nat_aqi_days['P45'][2100], positions = [xticks[2]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
# P37
p_37 =axsLeft[0, 0].boxplot(nat_aqi_days['P37'][2050], positions = [xticks[1] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))
axsLeft[0, 0].boxplot(nat_aqi_days['P37'][2100], positions = [xticks[2] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))

# Adjust parameters
axsLeft[0, 0].set_xlim([0, 8])
axsLeft[0, 0].set_xticks([])
axsLeft[0, 0].set_xticklabels([])
# axs[0, 0].set_xlabel('Year')
axsLeft[0, 0].set_ylabel('Days per year\n(U.S. national population-weighted average)', fontsize = lsize)
axsLeft[0, 0].set_title('Air Quality\n Advisories', fontsize = lsize)
axsLeft[0, 0].set_ylim([4, 25])
axsLeft[0, 0].set_yticks([5, 10, 15, 20, 25])
axsLeft[0, 0].set_yticklabels(['5', '10', '15', '20', '25'])
axsLeft[0, 0].spines['bottom'].set_visible(False)
axsLeft[0, 0].legend([r["boxes"][0], p_45["boxes"][0], p_37["boxes"][0]], ['REF', 'P45', 'P37'])

#### Create second plot
axsLeft[0, 1].boxplot(ad_2000, positions = [xticks[0]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
# AE = 1
axsLeft[0, 1].boxplot(nat_ad_days[1]['REF'][2050], positions = [xticks[1] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
axsLeft[0, 1].boxplot(nat_ad_days[1]['REF'][2100], positions = [xticks[2] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))

# P45
axsLeft[0, 1].boxplot(nat_ad_days[1]['P45'][2050], positions = [xticks[1]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
axsLeft[0, 1].boxplot(nat_ad_days[1]['P45'][2100], positions = [xticks[2]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
# P37
axsLeft[0, 1].boxplot(nat_ad_days[1]['P37'][2050], positions = [xticks[1] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))
axsLeft[0, 1].boxplot(nat_ad_days[1]['P37'][2100], positions = [xticks[2] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))

# Adjust parameters
axsLeft[0, 1].set_xlim([0, 8])
axsLeft[0, 1].set_xticks([])
# axs[0, 1].set_xticklabels(xticklabels)
# axs[0, 1].set_xlabel('Year')
axsLeft[0, 1].set_title('Optimal\n Adaptation', fontsize = lsize)
axsLeft[0, 1].set_ylim([198, None])
axsLeft[0, 1].spines['bottom'].set_visible(False)


##### Create 3rd
axsLeft[1, 0].boxplot(pm_high_2000, positions = [xticks[0]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
# REF
r = axsLeft[1, 0].boxplot(nat_aqi_days['REF'][2050], positions = [xticks[1] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
axsLeft[1, 0].boxplot(nat_aqi_days['REF'][2100], positions = [xticks[2] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))

# P45
p_45 = axsLeft[1, 0].boxplot(nat_aqi_days['P45'][2050], positions = [xticks[1]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
axsLeft[1, 0].boxplot(nat_aqi_days['P45'][2100], positions = [xticks[2]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
# P37
p_37 =axsLeft[1, 0].boxplot(nat_aqi_days['P37'][2050], positions = [xticks[1] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))
axsLeft[1, 0].boxplot(nat_aqi_days['P37'][2100], positions = [xticks[2] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))

# Adjust parameters
axsLeft[1, 0].set_xlim([0, 8])
axsLeft[1, 0].set_xticks(xticks)
axsLeft[1, 0].set_xticklabels(xticklabels)
axsLeft[1, 0].set_yticks([0])
axsLeft[1, 0].set_yticklabels(['0'])
axsLeft[1, 0].set_xlabel('Year', fontsize = lsize)
axsLeft[1, 0].set_ylim([0, 4])
axsLeft[1, 0].spines['top'].set_visible(False)
axsLeft[1, 0].tick_params(labeltop=False)

###### Create 4th plot
axsLeft[1, 1].boxplot(ad_2000, positions = [xticks[0]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
# AE = 1
axsLeft[1, 1].boxplot(nat_ad_days[1]['REF'][2050], positions = [xticks[1] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))
axsLeft[1, 1].boxplot(nat_ad_days[1]['REF'][2100], positions = [xticks[2] - width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[0], color= clist[0]),capprops=dict(color= clist[0]),
           whiskerprops=dict(color=clist[0]), medianprops=dict(color=clist[0]))

# P45
axsLeft[1, 1].boxplot(nat_ad_days[1]['P45'][2050], positions = [xticks[1]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
axsLeft[1, 1].boxplot(nat_ad_days[1]['P45'][2100], positions = [xticks[2]], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[1], color= clist[1]),capprops=dict(color= clist[1]),
           whiskerprops=dict(color=clist[1]), medianprops=dict(color=clist[1]))
# P37
axsLeft[1, 1].boxplot(nat_ad_days[1]['P37'][2050], positions = [xticks[1] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))
axsLeft[1, 1].boxplot(nat_ad_days[1]['P37'][2100], positions = [xticks[2] + width], widths = width, showfliers = False,
           patch_artist=True, boxprops=dict(facecolor= clist[2], color= clist[2]),capprops=dict(color= clist[2]),
           whiskerprops=dict(color=clist[2]), medianprops=dict(color=clist[2]))

# Adjust parameters
axsLeft[1, 1].set_xlim([0, 8])
axsLeft[1, 1].set_xticks(xticks)
axsLeft[1, 1].set_xticklabels(xticklabels)
axsLeft[1, 1].set_yticks([0])
axsLeft[1, 1].set_yticklabels(['0'])
axsLeft[1, 1].set_xlabel('Year', fontsize = lsize)
axsLeft[1, 1].set_ylim([0, 198])
axsLeft[1, 1].spines['top'].set_visible(False)
axsLeft[1, 1].tick_params(labeltop=False)
fig.subplots_adjust(hspace=0.05)


# Add cut lines
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
axsLeft[0, 0].plot([0, 1], [0, 0], transform=axsLeft[0, 0].transAxes, **kwargs)
axsLeft[1, 0].plot([0, 1], [1, 1], transform=axsLeft[1, 0].transAxes, **kwargs)
axsLeft[0, 1].plot([0, 1], [0, 0], transform=axsLeft[0, 1].transAxes, **kwargs)
axsLeft[1, 1].plot([0, 1], [1, 1], transform=axsLeft[1, 1].transAxes, **kwargs)

# Add a and b annotation
axsLeft[0,0].annotate('(a)', (6, 23.8), fontsize = lsize)
axsLeft[0,1].annotate('(b)', (6, 248), fontsize = lsize)





# Plot second subfigure
aqi_gdf_dict[2050].plot(color = aqi_gdf_dict[2050]['RGBA'], ax = axsMid[0])
aqi_gdf_dict[2100].plot(color = aqi_gdf_dict[2100]['RGBA'], ax = axsMid[1])

axsMid[0].set_title('P37 2050', fontsize = lsize)
axsMid[1].set_title('P37 2100', fontsize = lsize)

axsMid[0].annotate('(c)',(-125,27), fontsize = lsize)
axsMid[1].annotate('(d)',(-125,27), fontsize = lsize)

axsMid[0].set_axis_off()
axsMid[1].set_axis_off()



ax4 = subfigs[1].add_axes([0.05, .075, .9, 0.03])
cbar = mpl.colorbar.ColorbarBase(ax4, cmap = 'Oranges', orientation='horizontal', ticks = ticks, norm = norm)
cbar.set_label(label_aqi, fontsize = lsize)
cbar.ax.set_xticklabels(tick_labels)


#Plot third subfigurey
# Plot contour maps of adaptation days

gdf_dict['P37'][2050].plot(color = gdf_dict['P37'][2050]['RGBA'], ax = axsRight[0])
gdf_dict['P37'][2100].plot(color = gdf_dict['P37'][2100]['RGBA'], ax = axsRight[1])

axsRight[0].set_title('P37 2050', fontsize = lsize)
axsRight[1].set_title('P37 2100', fontsize = lsize)


axsRight[0].annotate('(e)',(-125,27), fontsize = lsize)
axsRight[1].annotate('(f)',(-125,27), fontsize = lsize)




axsRight[0].set_axis_off()
axsRight[1].set_axis_off()


ax5 = subfigs[2].add_axes([0.05, .075, .9, 0.03])
cbar = mpl.colorbar.ColorbarBase(ax5, cmap = 'Purples', orientation='horizontal', ticks = ticks, norm = norm)
cbar.set_label(label, fontsize = lsize)
cbar.ax.set_xticklabels(tick_labels)

title = 'fig2_combined_aqidays_adapt_contours_font'

plt.savefig(outdir + title + '.png', dpi = 400, bbox_inches = 'tight')






