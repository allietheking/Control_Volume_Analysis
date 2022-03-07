'''
Allie wrote this script back in 2020
It reads the definitions of the groups and their connectivity and plots them
In March 2022 added colors to indicate N/S/E/W
'''

########################################################################################
# import python packages
########################################################################################

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import geopandas as gpd
from shapely.ops import unary_union

#########################################################################################
# user input
#########################################################################################

# is it full resolution?
is_FR = True

# name for plots
if is_FR:
    plot_name = 'Full_Res'
else:
    plot_name = 'Agg'

# path to shapefile defining control volumes
if is_FR:
    cv_shapefile_path = r'X:\hpcshared\inputs\shapefiles\Agg_mod_contiguous_plus_subembayments_shoal_channel.shp'
else:
    cv_shapefile_path = r'X:\hpcshared\inputs\shapefiles\Agg_mod_contiguous_141.shp'

# path to control volumes definition file
if is_FR:
    control_volume_definition_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Control_Volume_Definitions\control_volume_definitions_FR.txt'
    connectivity_definitino_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Control_Volume_Definitions\connectivity_definitions_FR.txt'
else:
    control_volume_definition_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Control_Volume_Definitions\control_volume_definitions_141.txt'
    connectivity_definitino_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Control_Volume_Definitions\connectivity_definitions_141.txt'

# path for plots to check definitions
plot_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Plots\Group_Definitions'

# colors for N, S, E, W
colorN = 'r'
colorS = 'b'
colorW = 'm'
colorE = 'c'

###########################################################################################
# MAIN
###########################################################################################


# load dictionaries that define the groups (A, B, C, etc.) by the list of control volumes
# that make up each group (group_dict) and the list of pairs of control volumes making up
# fluxes out of the N, S, E, and W faces of each group (flux_pairs_dict)

# read from text file a dictionary grouping control volumes into lettered groups 
# format:
#   GROUP_NAME : CV1, CV2, CV3, ...
# where:
#   GROUP_NAME = name of group
#   CV1, CV2, CV3, ... = numbers of control volumes comprising the group
group_dict = {}
with open(control_volume_definition_path,'r') as f: 
    for line in f.readlines():
        # split line at ':' into key and control volume list, dropping newline character first
        (key, val) = line.rstrip('\n').split(':')  
        # remove whitespace surrounding key 
        key = key.strip()
        # add list to dictionary
        group_dict[key] = eval(val) 
      
        
# read from text file a dictionary specifying all the fluxes between 
# control volumes making up fluxes from each lettered group to the north, 
# south, east, and west
# format:
#   LETTER2DIRECTION : [[CVf,CVt]_1, [CVf,CVt]_2, [CVf,CVt]_3, ...]
# where:
#   LETTER = name of group
#   DIRECTION = N, S, E, or W stand for "to the north", south, east, or west, respectively
#   [CVf, CVt]_i = defines a flux from control volume CVf to control volume CVt
#   where i = 1, 2, 3, ... make up all the fluxes for a particular group/direction pair
flux_pairs_dict = {}
with open(connectivity_definitino_path,'r') as f: 
    for line in f.readlines():
        # split line at ':' into key and connectivity volume list, dropping newline character first
        (key, val) = line.rstrip('\n').split(':') 
        # trim whitespace around key
        key = key.strip()
        # add list to dictionary
        flux_pairs_dict[key] = eval(val) 

# load control volume shape file
gdf = gpd.read_file(cv_shapefile_path)

# loop through keys in group dictionary and plot the merged polygons 
for key in group_dict.keys():

    # this is hard coded ... for full res run, plot subset of control volumes relevant to the group at hand
    if is_FR:
        if '_Shoal_RMP_FR' in key or '_Channel_RMP_FR' in key:
            iplot = [168,169,170,155,157,161,158,159,160]
        elif '_RMP_FR' in key:
            iplot = [156,155,161,157,158,159,160]
        elif '_Shoal_WB_FR' in key or '_Channel_WB_FR' in key:
            iplot = [151,163,164,165,152,153,148,149,154,147]
        elif '_WB_FR' in key:
            iplot = [151,150,152,153,148,149,154,147]
        else: 
            iplot = np.arange(0,147)
    else:
        iplot = np.arange(0,len(gdf))

    # merge the control volume polygons comprising a group into a single polygon
    poly_list = []
    for cv in group_dict[key]:
        poly_list.append(gdf.iloc[cv]['geometry'])
    group_poly = unary_union(poly_list)
        
    # loop through the flux pairs dictionary and grab all the lists of flux pairs 
    # origninating in this group
    flux_pairs_keys = []
    for flux_pairs_key in flux_pairs_dict.keys():
        if flux_pairs_key[0:len(key)+3] == (key + ' to'):
            flux_pairs_keys.append(flux_pairs_key)
            
    # loop through all the lists of flux pairs originiating in this group and 
    # compile into one big list
    flux_pairs = []
    colors = []
    for flux_pairs_key in flux_pairs_keys:
        for flux_pairs_val in flux_pairs_dict[flux_pairs_key]:  
            # if a redundant flux pair is found, raise an error and exit, otherwise append
            # flux pair to the master list
            if flux_pairs_val in flux_pairs:
                print('Flux pair %s is redundant for group %s' % (flux_pairs_val, key))
                sys.exit()
            else:
                flux_pairs.append(flux_pairs_val)
                if 'to N' in flux_pairs_key:
                    colors.append(colorN)
                elif 'to S' in flux_pairs_key:
                    colors.append(colorS)
                elif 'to E' in flux_pairs_key:
                    colors.append(colorE)
                elif 'to W' in flux_pairs_key:
                    colors.append(colorW)
            
    # plot conktrol volume map with numbered control volumes
    fig = plt.figure(figsize=(40,30))
    ax = fig.subplots(1,1)
    gdf.iloc[iplot].plot(ax=ax,color='w',edgecolor='k')
    for cvnum in iplot:
        xc = gdf.iloc[cvnum].geometry.centroid.x 
        yc = gdf.iloc[cvnum].geometry.centroid.y
        plt.text(xc,yc,'%d' % cvnum, horizontalalignment='center',verticalalignment='center')
    
    # add the outline of the group polygon in red and put the title in red at the top
    ax.plot(*group_poly.exterior.xy,'r',linewidth=3)
    plt.title(key,color='r',fontsize=60)
    
    # loop through the flux pairs for this group, find the centroid of the starting 
    # and ending control volumes, and draw an error spanning middle 50% of centroid 
    # to centroid
    for flux_pair, color in zip(flux_pairs,colors):
    
        # starting and ending control volumes are defined by flux pair
        start_cv = flux_pair[0]
        end_cv = flux_pair[1]
        
        # find coordinates of centroids for starting and ending polygons
        xs = gdf.iloc[start_cv]['geometry'].centroid.x
        ys = gdf.iloc[start_cv]['geometry'].centroid.y
        xe = gdf.iloc[end_cv]['geometry'].centroid.x
        ye = gdf.iloc[end_cv]['geometry'].centroid.y
        
        # now plot an arrow representing the flux between control volumes
        plt.plot([xs, xe], [ys, ye], color=color, linewidth = 3)
        plt.plot(xe, ye, 'o', color=color, markersize=10)
        
    # save and close the figure
    plt.savefig(os.path.join(plot_path, plot_name + '_' + key + '.png'))
    plt.close()
    
        

    
