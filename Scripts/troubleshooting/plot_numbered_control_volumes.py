

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

# path to shapefile defining control volumes
cv_shapefile_path = '../Delta_Fullres_Polygons_Jabusch.shp'

###########################################################################################
# MAIN
###########################################################################################

# make directory for figures
if not os.path.exists('Check_Group_and_Connectivity_Definitions'):
    os.makedirs('Check_Group_and_Connectivity_Definitions')

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
with open('control_volume_definitions.txt','r') as f: 
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
with open('connectivity_definitions.txt','r') as f: 
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
    for flux_pairs_key in flux_pairs_keys:
        for flux_pairs_val in flux_pairs_dict[flux_pairs_key]:  
            # if a redundant flux pair is found, raise an error and exit, otherwise append
            # flux pair to the master list
            if flux_pairs_val in flux_pairs:
                print('Flux pair %s is redundant for group %s' % (flux_pairs_val, key))
                sys.exit()
            else:
                flux_pairs.append(flux_pairs_val)
            
    # plot conktrol volume map with numbered control volumes
    fig = plt.figure(figsize=(40,30))
    ax = fig.subplots(1,1)
    gdf.plot(ax=ax,color='w',edgecolor='k')
    for cvnum in range(len(gdf)):
        xc = gdf.iloc[cvnum].geometry.centroid.x 
        yc = gdf.iloc[cvnum].geometry.centroid.y
        plt.text(xc,yc,'%d' % cvnum, horizontalalignment='center',verticalalignment='center')
    
    # add the outline of the group polygon in red and put the title in red at the top
    ax.plot(*group_poly.exterior.xy,'r',linewidth=3)
    plt.title(key,color='r',fontsize=60)
    
    # loop through the flux pairs for this group, find the centroid of the starting 
    # and ending control volumes, and draw an error spanning middle 50% of centroid 
    # to centroid
    for flux_pair in flux_pairs:
    
        # starting and ending control volumes are defined by flux pair
        start_cv = flux_pair[0]
        end_cv = flux_pair[1]
        
        # find coordinates of centroids for starting and ending polygons
        xs = gdf.iloc[start_cv]['geometry'].centroid.x
        ys = gdf.iloc[start_cv]['geometry'].centroid.y
        xe = gdf.iloc[end_cv]['geometry'].centroid.x
        ye = gdf.iloc[end_cv]['geometry'].centroid.y
        
        # now plot an arrow representing the flux between control volumes
        plt.plot([xs, xe], [ys, ye], 'm', linewidth = 3)
        plt.plot(xe, ye, 'mo', markersize=10)
        
    # save and close the figure
    plt.savefig('Check_Group_and_Connectivity_Definitions/' + key + '.png')
    plt.close()
    
        

    
