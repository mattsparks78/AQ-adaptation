# -*- coding: utf-8 -*-
"""
Created on Fri May 13 10:50:05 2022

@author: User
"""


import pandas as pd
import pickle as cp


# Load benmap grid
AQgrid = pd.read_csv(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\BenMapGrid.csv')
AQgrid.drop(218, inplace = True)

# Create list of row,county pairs in grid cell
rc = list(zip(list(AQgrid['ROW']),list(AQgrid['COL'])))

# Load Counties to AQgrid
dtype_dict = {'ROW':int, 'COL': int, 'FIPS': str, 'pop_share':float}
counties_grid = pd.read_csv(r'C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\county_population_per_cell_b.csv', dtype = dtype_dict)

# Create dictionary to hold population share data
pop_share = {}
for r, c in rc:
    pop_share[(r, c)] = {}
    gc = counties_grid.loc[(counties_grid['ROW'] == r)*(counties_grid['COL'] == c)]
    for j in gc['FIPS']:
        pop_share[(r, c)][j] = gc.loc[gc['FIPS'] == j, 'pop_share'].item()
        
# Create dictionary to hold SCI values
sci_dict = {}
for i in rc:
    sci_dict[i] = {}
    for j in rc:
        sci_dict[i][j] = 0

# Load FB county to county SCI data
# ctc = pd.read_csv(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\SCI\SCI_county_county.tsv", dtype = 'object', sep = '\t')

# Find unique counties in ctc to create indices in dictionary
# u_c = ctc.user_loc.unique()

# Create dictionary for county to county SCI data
# ctc_dict = {}
# for i in u_c:
#     ctc_dict[str(i)] = {}

# # Add SCI data from county to county data frame to dict
# for i in u_c:
#     ctc_dict[str(i)] = {}
#     cncts = ctc.loc[ctc['user_loc'] == i] 
#     for j in cncts['fr_loc']:
#         ctc_dict[str(i)][str(j)] = int(cncts.loc[cncts['fr_loc'] == j, 'scaled_sci'].item())

# #Output ctc dict to pickle        
# with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\ctc_sci", 'wb') as f:
#      cp.dump(ctc_dict,f)
# f.close()

# Import ctc dict (previously calculated with code above)
with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\BenMap Grid\ctc_sci", 'rb') as f:
     ctc_dict = cp.load(f)
f.close()        

i = 0
# Loop over all grid cells
for r,c in rc:
    # Isolate current user grid cell 
    uc = pop_share[(r, c)]
        
    # Loop over all friend cells. 
    for row, col in rc:
        # Isolate current friend cell
        fc = pop_share[(row, col)]
        # Create list for user-friend sci values for the grid cell. Will add up pop-weighted county-county values
        fr_sci = []
        
        # Iterate over counties in user cell
        for fips in uc:
            # Select list of SCI values based on user county
            if fips not in ['51515', '51560', '46113']:
                c_sci = ctc_dict[fips]

            # Iterate over counties in friend cell
                for fr_county in fc:
                   if fr_county in c_sci:
                       scen_sci = c_sci[fr_county]
                       pop_share_base = uc[fips]
                       pop_share_friend = fc[fr_county]
                       fr_sci.append(round(scen_sci * pop_share_base * pop_share_friend))
                   else:
                       fr_sci.append(1)
                  
        sci_dict[(r, c)][(row, col)] = sum(fr_sci)
    i = i+1
    print(i)

# save sci_dict
with open(r"C:\Users\User\OneDrive - University of Waterloo\Sparks\PNAS\WorkingFolder\Data\SCI\sci_tot_b", 'wb') as f:
      cp.dump(sci_dict,f)
f.close()














                
            