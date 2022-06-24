
# make some plots to check for mass conservation, using the results of step3
# note this checks the reactions only, for composite parameters
# Allie April 2022
#

#################################################################
# IMPORT PACKAGES
#################################################################

import sys, os, re
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
import logging
import datetime
import geopandas as gpd
import socket
hostname = socket.gethostname()
if hostname == 'richmond':
    raise Exception("Do not run this script on richmond until we update the conda environment... run on chicago or your laptop instead")


#################################################################
# USER INPUT
#################################################################

# which parameter to check?
#param = 'TN'
param = 'TN_include_sediment'

# which reaction rate to normalize by?
normalize_by = 'NO3,dDenit'

# run id 
#runid = 'FR16_28'
#runid = 'FR13_003'
#runid = 'FR13_023'
#runid = 'FR17_003'
#runid = 'FR17_017'
runid = 'FR18_005'

# is it the delta? if so, assumes full resolution
is_delta = False

# path to base level balance tables (organized by runid within this folder)
balance_table_path = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables'

# output plot path
if is_delta:
    balance_table_id = 'Delta_' + runid
else:
    balance_table_id = runid
output_path = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Plots/%s/mass_conservation_check/' % balance_table_id

# name of the pdf file that will contain maps showing contribution of the grouped reaction terms
rx_map_path = os.path.join(output_path, '%s_%s_reaction_percentage_map.pdf' % (runid,param))

# root directory (for runs, shapefiles, and output csv files)
root_dir = '/richmondvol1/hpcshared/'   ## running on chicago or richmond
#root_dir = r'X:\hpcshared'              ## running on windows laptop with mounted drive

# is full resolution?
if 'FR' in runid:
    FR = True
else:
    FR = False

# path to shapefiles corresponding to the polygons in the balance tables, also select a subset of polygons to plot
if is_delta:
    shpfn_poly = os.path.join(root_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Polygons_Dave_Plus_WB_v4.shp')  
    rx_pos = (508800,4286000) 
    balance_table_id = 'Delta_' + runid

elif FR:
    shpfn_poly = os.path.join(root_dir,'inputs','shapefiles','Agg_mod_contiguous_plus_subembayments_shoal_channel.shp')
    iplot = [0, 117, 139, 2, 113, 114, 115, 111, 1, 3, 116, 7, 4, 112, 5, 108, 109, 110, 9, 107, 
    6, 8, 29, 10, 12, 11, 138, 137, 100, 19, 18, 13, 20, 26, 25, 21, 14, 24, 101, 102, 103, 104, 
    105, 106, 28, 27, 22, 15, 30, 23, 35, 34, 33, 32, 31, 17, 41, 40, 39, 38, 37, 36, 16, 140, 
    44, 43, 42, 46, 45, 49, 47, 97, 98, 86, 96, 95, 85, 93, 52, 87, 94, 48, 51, 50, 90, 88, 91, 
    89, 92, 53, 56, 99, 54, 55, 57, 67, 65, 66, 136, 58, 59, 63, 64, 60, 62, 61, 143, 144, 141, 
    68, 69, 70, 71, 78, 72, 84, 79, 146, 80, 82, 83, 142, 145, 73, 81, 74, 76, 77, 75]
    rx_pos = (573900,4204500)
    balance_table_id = runid
else:
    shpfn_poly =  os.path.join(root_dir,'inputs','shapefiles','Agg_mod_contiguous_141.shp')
    iplot = [  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
        13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
        26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
        39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
        52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
        65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
        78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
        91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
       134, 136, 137, 138, 139, 140]
    rx_pos = (573900,4204500)
    balance_table_id = runid

#################################################################
# FUNCTIONS
#################################################################

# define a function that returns a list of the reactions included in a balance table
def get_rx_list(df, substance):    
    rx_list = []
    columns = df.columns  
    for column in columns:
        if not column in ['time', 'Control Volume', 'Concentration (mg/l)', 'Volume',
                          'Volume (Mean)', 'Area', 'dVar/dt', '%s,Loads in' % substance, '%s,Loads out' % substance,
                          '%s,Transp in' % substance, '%s,Transp out' % substance]:
            if not 'Flux' in column:
                if not 'To_poly' in column:
                    rx_list.append(column)
    return rx_list


#################################################################
# MAIN
#################################################################

# make directory for output plots
if not os.path.exists(output_path):
    os.makedirs(output_path)

# read table
table_name = '%s_Table.csv' % param.lower()
df = pd.read_csv(os.path.join(balance_table_path,balance_table_id,table_name))

# read the polygon shapefile, and if delta, find the water board polygons and select those to plot
gdf = gpd.read_file(shpfn_poly)
if is_delta:
    ind = gdf['polytype'] == 'WB'
    iplot = list(gdf.loc[ind].index)

# select subset of polygons to plot, and come up with list of corresponding polygon names
iplot = np.sort(iplot)
gdf = gdf.iloc[iplot]
polys = np.array(['polygon%d' % i for i in iplot])

# get list of reactions
rx_list = get_rx_list(df, param)

# plot style list
plot_colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000',
               '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000',
               '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
line_styles = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
               '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--', '--',
               ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':', ':']

# initialize a dataframee with reactions as columns and polygons as rows
df_percent = pd.DataFrame(columns=polys, index=pd.Index(rx_list))

# loop through polygons, plot reaction terms, and tabulate errors for terms that should be zero
for poly in polys:

    # get data just for this polygon
    ind = df['Control Volume']==poly
    df1 = df.loc[ind]
    time = df1['time'].astype('datetime64[ns]')

    # create figure with all reactions
    fig, ax = plt.subplots(figsize=(16,16))
    counter2 = 0
    for rx in rx_list:
        ax.plot(time, df1[rx], label=rx, color = plot_colors[counter2], linestyle = line_styles[counter2])
        counter2 = counter2 + 1
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
    ax.set_ylabel('reaction rate g/d')
    ax.set_title('%s: %s' % (runid, poly))
    fig.tight_layout()
    plt.savefig(os.path.join(output_path,'%s_%s_%s_check_rx_mass_cons.png' % (runid, param, poly)))
    plt.close('all')

    # compute time average of all reaction terms
    df1_tavg = df1[rx_list].mean(axis=0)

    # compute the percent of each reaction compared to the absolute value of the "normalize by" reaction
    df1_percent = df1_tavg/np.abs(df1_tavg[normalize_by])*100
    
    # tabulate percent
    df_percent[poly] = df1_percent

# now loop through the reaction terms and create a pdf that shows percent contribution of each, w.r.t. max souce/sink
with PdfPages(rx_map_path) as pdf:

    for rx in rx_list:

        # load up the geodataframe with the percent values correcponding to this reaction
        gdf['percent'] = df_percent.loc[rx].values

        # format the reaction name for printing on the plot
        rx1 = rx.replace(' + ',' +\n              ')

        # plot percent
        fig, ax = plt.subplots(figsize=(11, 11))
        gdf.plot(ax=ax,column='percent',vmin=-100,vmax=100,cmap='seismic',legend=True)
        gdf.boundary.plot(ax=ax,color='k')
        ax.axis('off')
        ax.set_title('%s: reaction rate as percent of |%s|\ntime averaged over entire simulation' % (runid, normalize_by))
        ax.text(rx_pos[0],rx_pos[1],rx1,ha='left',va='top',fontsize=14)
        fig.tight_layout()
        pdf.savefig()
        plt.close('all')



