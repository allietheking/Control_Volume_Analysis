

########################################################################################
# import python packages
########################################################################################

import sys
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import geopandas as gpd
from shapely.ops import unary_union

#########################################################################################
# user input
#########################################################################################

# path to shapefile defining control volumes
cv_shapefile_path = 'Agg_mod_contiguous_141.shp'

# clusters of groups for making maps
group_cluster_dict = {'Groups' : ['Y','Z','W','X','S1','S2','T2','T1','U2','U1','V','A','B','C','D','E','F','G','H','I','J','K','L'],
                      'RMP Subembayments' : ['SB_RMP_west_shoal', 'SB_RMP_channel', 'SB_RMP_east_shoal', 'LSB'],
                      'WBR2 Subembayments' : ['SB_WB_west_shoal', 'SB_WB_channel', 'SB_WB_east_shoal', 'LSB']}
                      
# group label dictionary
group_label_dict = {'Y' : 'Y','Z' : 'Z','W' : 'W','X' : 'X',
                   'S1' : 'S1','S2' : 'S2','T2' : 'T2','T1' : 'T1','U2' : 'U2','U1' : 'U1',
                   'V' : 'V','A' : 'A','B' : 'B','C' : 'C','D' : 'D',
                   'E' : 'E','F' : 'F','G' : 'G','H' : 'H','I' : 'I',
                   'J' : 'J','K' : 'K','L' : 'L', 
                   'LSB' : 'Lower\nSouth\nBay',
                   'SB_RMP_west_shoal' : 'South Bay\nWest Shoal',
                   'SB_RMP_channel' : 'South Bay\nChannel',
                   'SB_RMP_east_shoal' : 'South Bay\nEast Shoal',
                   'SB_WB_west_shoal' : 'South Bay\nWest Shoal',
                   'SB_WB_channel' : 'South Bay\nChannel',
                   'SB_WB_east_shoal' : 'South Bay\nEast Shoal'}
                   
# group label coordinate dictionary
group_label_coord_dict = {'LSB' : (580696.9672299023, 4147935.943827709),
                   'SB_RMP_west_shoal' : (561492.347787714, 4153519.6657800633),
                   #'SB_RMP_channel' : 'South Bay\nChannel',
                   'SB_RMP_east_shoal' : (574013.4212566298, 4165871.535553453),
                   'SB_WB_west_shoal' : (553624.3759457602, 4160795.4246876766),
                   #'SB_WB_channel' : 'South Bay\nChannel',
                   'SB_WB_east_shoal' : (572998.1990834743, 4169424.813159497)}


# zoom window for maps
zoom_window = (547227.2676570563, 593439.214100827, 4136961.3921359004, 4198717.356928939)

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
with open('control_volume_definitions_141.txt','r') as f: 
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
with open('connectivity_definitions_141.txt','r') as f: 
    for line in f.readlines():
        # split line at ':' into key and connectivity volume list, dropping newline character first
        (key, val) = line.rstrip('\n').split(':') 
        # trim whitespace around key
        key = key.strip()
        # add list to dictionary
        flux_pairs_dict[key] = eval(val) 

# load control volume shape file
gdf = gpd.read_file(cv_shapefile_path)

# plot original control volumes
fig = plt.figure(figsize=(10,8))
ax = fig.subplots(1,1)
gdf.plot(ax=ax,color='gray')
for cv in range(len(gdf)):
    xc = gdf.iloc[cv].geometry.centroid.x
    yc = gdf.iloc[cv].geometry.centroid.y
    plt.text(xc,yc,'%d' % cv, horizontalalignment='center',verticalalignment='center', color='k', fontsize=6)
plt.savefig('CV_Map.png')

# initialize list of groups and their geometries
groups_list = []
geometry_list = []

# loop through keys in group dictionary and plot the merged polygons 
for key in group_dict.keys():

    # merge the control volume polygons comprising a group into a single polygon
    poly_list = []
    for cv in group_dict[key]:
        poly_list.append(gdf.iloc[cv]['geometry'])
    group_poly = unary_union(poly_list)
    
    # append lists of groups and their polygons
    groups_list.append(key)
    geometry_list.append(group_poly)
    
# create geopandas dataframe for groups
df = pd.DataFrame(index=groups_list)
gdf_groups = gpd.GeoDataFrame(df, geometry=geometry_list)

# loop through group clusters
for igc in range(len(group_cluster_dict)):

    # name of gropup cluster
    group_cluster_name = list(group_cluster_dict.keys())[igc]

    # select subset of dataframe with group cluster groups
    gdf_cluster = gdf_groups.loc[group_cluster_dict[group_cluster_name]]    

    # plot control volume map with numbered control volumes
    fig = plt.figure(figsize=(10,12))
    ax = fig.subplots(1,1)
    gdf.plot(ax=ax,color='gray')
    gdf_cluster.geometry.boundary.plot(ax=ax,edgecolor='k')
    for group in gdf_cluster.index:
        if group in group_label_coord_dict.keys():
            xc = group_label_coord_dict[group][0]
            yc = group_label_coord_dict[group][1]
        else:
            xc = gdf_cluster.loc[group].geometry.centroid.x
            yc = gdf_cluster.loc[group].geometry.centroid.y
        plt.text(xc,yc,'%s' % group_label_dict[group], horizontalalignment='center',verticalalignment='center', color='k', fontsize=16)
    plt.title(group_cluster_name, fontsize=20)
    plt.axis(zoom_window)
    plt.savefig('Group_Map_' + group_cluster_name.replace(' ','_') + '.png')
        
   

    
