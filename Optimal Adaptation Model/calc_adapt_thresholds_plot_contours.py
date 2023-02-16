# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:09:50 2022

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
indir1 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\General\\"
indir2 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Future AQ\\"
indir3 = r"C:\Users\User\OneDrive - University of Waterloo\ClimatePenalty\WorkingFolder\Data\AGU2020-Fig3\\"
indir4 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\US Boundary\\"
indir5 = r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\cb_2018_us_nation_5m\\'
indir6 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\\"
indir7 = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\Rational Actor\\"
outdir = r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Figures\\"

# Load future concentration data. REF, P3.7 or P4.5
with open(indir2 + 'REF_PM25_2050_2100_d', 'rb') as f:
    pm_out = cp.load(f)
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


# Import baseline mortality per grid cell
bm = pd.read_csv(indir1 + 'baseline_mortality_b.csv')

# Create list of all grid cells
rc = list(zip(list(wage['ROW']), list(wage['COL'])))

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

# Define US boundary
us_bound = gpd.read_file(indir5 + 'cb_2018_us_nation_5m.shp')
us_bound = us_bound.to_crs(epsg=4326)

# Calculate mean pm concentrations for each grid cell and threshold for adaptation
pm_cell_mean = {}
thresh_cell = {}

for year in [2050, 2100]:
    pm_cell_mean[year] = {}
    thresh_cell[year] = {}
    if year == 2050:
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
    
    for pol in ['REF', 'P45', 'P37']:
        pm_cell_mean[year][pol] = {}
        thresh_cell[year][pol] = {}
        fname = 'pm_base_' + pol + '_cell_' + str(year)
        with open(indir7 + fname, 'rb') as f:
             pm_cell = cp.load(f)
        f.close()
        for r, c in rc:
            pm_cell_mean[year][pol][(r, c)] = np.mean(pm_cell[(r, c)])
            w_i = wage.loc[(wage['ROW'] == r)&(wage['COL'] == c), '2000'].item() * EPF
            base_mort = bm.loc[(bm['ROW'] == r)&(bm['COL'] == c), 'BaseMortRate'].item() * PMMF
            dPM = (AT*w_i*10)/(base_mort*AF*VSL)
            thresh_cell[year][pol][(r, c)] = 365 * pm_cell_mean[year][pol][(r, c)] - 364 * (pm_cell_mean[year][pol][(r, c)] - dPM)

# Convert PM threshold to AQI
thresh_cell_AQI = {}
for year in [2050, 2100]:
    thresh_cell_AQI[year] = {}
    for pol in ['REF', 'P45', 'P37']:
        thresh_cell_AQI[year][pol] = {}
        for r, c in rc:
            con = thresh_cell[year][pol][(r, c)]
            if con <= 12:
                aqi = round(4.17*con)
            elif 12 < con <= 35.4:
                aqi = round(2.1*(con - 12) + 51)
            elif 35.4 < con <= 55.4:
                aqi = round(2.46*(con - 35.4) + 101)
            elif 55.4 < con <= 150.4:
                aqi = round(0.52*(con-55.4) + 151)
            elif 150.4 < con <= 250.4:
                aqi = round(0.99*(con-150.4) + 201)
            else:
                aqi = round(0.8*(con-250.4) + 301)
            thresh_cell_AQI[year][pol][(r, c)] = aqi

#%%
#Plot contours of AQI thresholds
grid = gpd.read_file(indir6 + "BenMapGridPoints.csv")
grid = grid.drop('geometry', axis = 1)
grid['geom'] = grid['geom'].apply(wkt.loads)
grid = grid.set_geometry('geom')
grid = grid.set_crs(epsg=4326)
grid['ROW'] = grid['ROW'].astype('int')
grid['COL'] = grid['COL'].astype('int')
grid = grid.drop([1, 219, 220])

# Plot adaptation days by grid cell
aqi_gdf = copy.deepcopy(grid)
for pol in ['P37', 'P45', 'REF']:
    for year in [2050, 2100]:
        key = pol + ' ' + str(year)
        aqi_gdf[key] = 0

for pol in ['P37', 'P45', 'REF']:
    for year in [2050, 2100]:
        key = pol + ' ' + str(year)
        if year == 2050:
            for r, c in rc:
                aqi_gdf.loc[(aqi_gdf['ROW'] == r) & (aqi_gdf['COL'] == c), key] = thresh_cell_AQI[year][pol][(r, c)]
        else:
            for r, c in rc:
                aqi_gdf.loc[(aqi_gdf['ROW'] == r) & (aqi_gdf['COL'] == c), key] = thresh_cell_AQI[year][pol][(r, c)] 


mymax = np.max([np.max(aqi_gdf['P45 2100']),np.max(aqi_gdf['P37 2100']),np.max(aqi_gdf['REF 2100'])])

norm = mpl.colors.Normalize(vmin= 0, vmax= mymax)

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

gdf_dict = {}
for pol in ['P37', 'P45', 'REF']:
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
                b = aqi_gdf.loc[(aqi_gdf['Lon'] == str(lon_arr[j])) & (aqi_gdf['Lat'] == str(lat_arr[i])), key]
                if len(b) == 0:
                    pass 
                else:
                    ben_arr[i, j] = aqi_gdf.loc[(aqi_gdf['Lon'] == str(lon_arr[j])) & (aqi_gdf['Lat'] == str(lat_arr[i])), key].item() 
        
        lon_arr = sp.ndimage.zoom(lon_arr, 5)
        lat_arr = sp.ndimage.zoom(lat_arr, 5)
        ben_arr = sp.ndimage.zoom(ben_arr, 5)
        
                
        collec_poly = plt.contourf(lon_arr, lat_arr, ben_arr, levels = 25, cmap = 'Reds_r', norm = norm);
        
        test = collec_to_gdf(collec_poly)
        
        ticks=[0, mymax*1/6,mymax*1/3,mymax*1/2,mymax*2/3,mymax*5/6,mymax]
        label = 'AQI threshold at which single-day optimal adaptation yields net benefits'
        tick_labels = [str(0), str(round(mymax*1/6)),str(round(mymax/3)),str(round(mymax/2)),str(round(mymax*2/3)),str(round(mymax*5/6)),str(round(mymax))]
        
        
        colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
        test['RGBA'] = colors
        test = test.to_crs(epsg=4326)
        for i in range(len(test.geometry)):
            test['geometry'][i] = test['geometry'][i].intersection(us_bound['geometry'][0])
        
        gdf_dict[pol][year] = test

#%% Plot the contour maps
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6))= plt.subplots(2, 3, figsize = (10, 6))

gdf_dict['P37'][2050].plot(color = gdf_dict['P37'][2050]['RGBA'], ax = ax1)
gdf_dict['P45'][2050].plot(color = gdf_dict['P45'][2050]['RGBA'], ax = ax2)
gdf_dict['REF'][2050].plot(color = gdf_dict['REF'][2050]['RGBA'], ax = ax3)

gdf_dict['P37'][2100].plot(color = gdf_dict['P37'][2100]['RGBA'], ax = ax4)
gdf_dict['P45'][2100].plot(color = gdf_dict['P45'][2100]['RGBA'], ax = ax5)
gdf_dict['REF'][2100].plot(color = gdf_dict['REF'][2100]['RGBA'], ax = ax6)

# fig.suptitle('AQI Threshold for Adaptation', fontsize=16)

ax1.set_title('P37 2050', fontsize = 10)
ax2.set_title('P45 2050', fontsize = 10)
ax3.set_title('REF 2050', fontsize = 10)

ax4.set_title('P37 2100', fontsize = 10)
ax5.set_title('P45 2100', fontsize = 10)
ax6.set_title('REF 2100', fontsize = 10)

ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()
ax4.set_axis_off()
ax5.set_axis_off()
ax6.set_axis_off()

ax7 = fig.add_axes([0.2, .075, .6, 0.03])
cbar = mpl.colorbar.ColorbarBase(ax7, cmap = 'Reds_r', orientation='horizontal', ticks = ticks, norm = norm)
cbar.set_label(label)
cbar.ax.set_xticklabels(tick_labels)

fig.tight_layout()
plt.show()

fig.savefig(outdir + 'AQI thresholds all scen contour.png', dpi = 400, bbox_inches = 'tight')
