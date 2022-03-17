
''' allie'''

##################
# IMPORT MODULES
##################

import copy
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
import cmocean 
try:
    import geopandas as gpd
except:
    raise Exception('\n' + 
          'ERROR: this Python environment does not include geopandas.\n' + 
          'To load environment with geopandas, exit Python, execute the\n' + 
          'following at the command line:\n' + 
          '    source activate /home/zhenlin/.conda/envs/my_root\n' + 
          'and then re-launch Python.\n')
    

###################
## USER INPUT
###################


# figure size
figsize = (5,5)

# run ID
#run_ID = 'FR13_003'
run_ID = 'FR17_003'
#run_ID = 'FR13_021'
#run_ID = 'FR17_014'
#run_ID = 'FR18_004'

# water year
#wy = 2013
wy = 2017
#wy = 2018

# what rate to plot? options are "denit", 'dpp', "nloss", "oxycon"
rate_name = 'dpp'

# time period (Annual, Seasonal, Monthly, maybe add Weekly later?)
time_period = 'Monthly'

flag = 5
if flag==0:
    rate_name = 'denit'
    time_period = 'Monthly'
elif flag==1:
    rate_name = 'denit'
    time_period = 'Seasonal'
elif flag==2:
    rate_name = 'dpp'
    time_period = 'Monthly'
elif flag==3:
    rate_name = 'dpp'
    time_period = 'Seasonal'
elif flag==4:
    rate_name = 'nloss'
    time_period = 'Monthly'
elif flag==5:
    rate_name = 'nloss'
    time_period = 'Seasonal'

# how to norm
norm = 'Area'  # options are 'None', 'Area', 'Volume'

# path to balance tables (they should be organized by run ID within this folderr)
balance_table_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Balance_Tables\%s' % run_ID

# output folder path
figure_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Plots\%s\rate_maps' % run_ID

# path to the shapefile w/ the base level control volumes
if 'FR' in run_ID:
    shp_fn = r'X:\hpcshared\inputs\shapefiles\Agg_mod_contiguous_plus_subembayments_shoal_channel.shp'
    iplot = [0, 117, 139, 2, 113, 114, 115, 111, 1, 3, 116, 7, 4, 112, 5, 108, 109, 110, 9, 107, 
    6, 8, 29, 10, 12, 11, 138, 137, 100, 19, 18, 13, 20, 26, 25, 21, 14, 24, 101, 102, 103, 104, 
    105, 106, 28, 27, 22, 15, 30, 23, 35, 34, 33, 32, 31, 17, 41, 40, 39, 38, 37, 36, 16, 140, 
    44, 43, 42, 46, 45, 49, 47, 97, 98, 86, 96, 95, 85, 93, 52, 87, 94, 48, 51, 50, 90, 88, 91, 
    89, 92, 53, 56, 99, 54, 55, 57, 67, 65, 66, 136, 58, 59, 63, 64, 60, 62, 61, 143, 144, 141, 
    68, 69, 70, 71, 78, 72, 84, 79, 146, 80, 82, 83, 142, 145, 73, 81, 74, 76, 77, 75]
else:
    shp_fn = r'X:\hpcshared\inputs\shapefiles\Agg_mod_contiguous_141.shp'
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

# figure output path and name
figure_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Plots\%s\rate_maps' % run_ID

# nanpercentile for color map cutoff
cper = 97.5

# here's where we set up what goes into these different rates
if rate_name=='denit':
    
    # title for figure
    rate_title = 'Denitrification'

    # grams of what in the units?   
    grams_of_what = 'N'

    # list of balance tables
    balance_table_list = ['no3_Table.csv']

    # multiplier for each balance table
    multiplier_list = [-1]

    # for each balance table, list of reactions to sum
    reaction_list = [['NO3,dDenitSed', 'NO3,dNiDen']]

    # color maps
    cmap = 'Spectral_r'

    # center at zero?
    cmap_diverging = False

    # include a log scale plot?
    log_scale = True

elif rate_name=='dpp':

    # title for figure
    rate_title = 'Net Primary Productivity'

    # grams of what in the units?   
    grams_of_what = 'C'

    # old version of the model doesn't have greens 
    if run_ID in ['FR13_003','FR17_003']:

        # list of balance tables
        balance_table_list = ['diat_Table.csv']

        # multiplier for each balance table
        multiplier_list = [1]

        # for each balance table, list of reactions to sum
        reaction_list = [['Diat,dPPDiat','Diat,dcPPDiat']]

    # newer versions do have greens
    else:

        # list of balance tables
        balance_table_list = ['diat_Table.csv', 'green_Table.csv']

        # multiplier for each balance table
        multiplier_list = [1, 1]

        # for each balance table, list of reactions to sum
        reaction_list = [['Diat,dPPDiat','Diat,dcPPDiat'],
                         ['Green,dPPGreen','Green,dcPPGreen']]
    # color maps
    cmap = cmocean.cm.algae

    # center at zero?
    cmap_diverging = False

    # include a log scale plot?
    log_scale = True

elif rate_name=='nloss':

    # title for figure
    rate_title = 'TN Reactive Loss'

    # grams of what in the units?   
    grams_of_what = 'N'

    # old version of the model doesn't have greens and diatom multiplier is 0.15
    if run_ID in ['FR13_003','FR17_003']:

        # list of balance tables
        balance_table_list = ['no3_Table.csv', 'nh4_Table.csv', 'pon1_Table.csv', 'don_Table.csv', 'diat_Table.csv', 'zoopl_v_Table.csv', 'zoopl_r_Table.csv', 'zoopl_e_Table.csv']

        # multiplier for each balance table
        multiplier_list = [-1, -1, -1, -1, -0.15, -0.1818, -0.1818, -0.1818]

        # for each balance table, list of reactions to sum
        reaction_list = [['NO3,dDenitWat', 'NO3,dNitrif', 'NO3,dDenitSed', 'NO3,dNiDen', 'NO3,dNO3Upt'],
                     ['NH4,dMinPON1', 'NH4,dMinDON', 'NH4,dNitrif', 'NH4,dMinDetNS1',
                      'NH4,dMinDetNS2', 'NH4,dZ_NRes', 'NH4,dNH4Aut', 'NH4,dNH4Upt'], 
                     ['PON1,dCnvPPON1', 'PON1,dCnvDPON1', 'PON1,dMinPON1', 'PON1,dZ_NMrt',
                      'PON1,dZ_NDef', 'PON1,dZ_NSpDet', 'PON1,dZ_PON1', 'PON1,dSedPON1',
                      'PON1,dMortDetN', 'PON1,dResS1DetN', 'PON1,dResS2DetN',
                      'PON1,dResS1DiDN', 'PON1,dResS2DiDN'], 
                     ['DON,dCnvDPON1', 'DON,dMinDON'],
                     ['Diat,dZ_Diat', 'Diat,dSedDiat', 'Diat,dPPDiat', 'Diat,dMrtDiat','Diat,dcPPDiat'],
                     ['Zoopl_V,dZ_Vgr', 'Zoopl_V,dZ_Vmor'],
                     ['Zoopl_R,dZ_SpwDet', 'Zoopl_R,dZ_Rgr', 'Zoopl_R,dZ_Rmor'],
                     ['Zoopl_E,dZ_Ea', 'Zoopl_E,dZ_Ec', 'Zoopl_E,dZ_Emor']]

    # newer versions of the model do have greens and diatom multiplier is 0.16
    else: 

        # list of balance tables
        balance_table_list = ['no3_Table.csv', 'nh4_Table.csv', 'pon1_Table.csv', 'don_Table.csv', 'diat_Table.csv', 'green_Table.csv', 'zoopl_v_Table.csv', 'zoopl_r_Table.csv', 'zoopl_e_Table.csv']

        # multiplier for each balance table
        multiplier_list = [-1, -1, -1, -1, -0.16, -0.16, -0.1818, -0.1818, -0.1818]

        # for each balance table, list of reactions to sum
        reaction_list = [['NO3,dDenitWat', 'NO3,dNitrif', 'NO3,dDenitSed', 'NO3,dNiDen', 'NO3,dNO3Upt'],
                     ['NH4,dMinPON1', 'NH4,dMinDON', 'NH4,dNitrif', 'NH4,dMinDetNS1',
                      'NH4,dMinDetNS2', 'NH4,dZ_NRes', 'NH4,dNH4Aut', 'NH4,dNH4Upt'], 
                     ['PON1,dCnvPPON1', 'PON1,dCnvDPON1', 'PON1,dMinPON1', 'PON1,dZ_NMrt',
                      'PON1,dZ_NDef', 'PON1,dZ_NSpDet', 'PON1,dZ_PON1', 'PON1,dSedPON1',
                      'PON1,dMortDetN', 'PON1,dResS1DetN', 'PON1,dResS2DetN',
                      'PON1,dResS1DiDN', 'PON1,dResS2DiDN'], 
                     ['DON,dCnvDPON1', 'DON,dMinDON'],
                     ['Diat,dZ_Diat', 'Diat,dSedDiat', 'Diat,dPPDiat', 'Diat,dMrtDiat','Diat,dcPPDiat'],
                     ['Green,dZ_Grn', 'Green,dSedGreen', 'Green,dPPGreen', 'Green,dMrtGreen','Green,dcPPGreen'],
                     ['Zoopl_V,dZ_Vgr', 'Zoopl_V,dZ_Vmor'],
                     ['Zoopl_R,dZ_SpwDet', 'Zoopl_R,dZ_Rgr', 'Zoopl_R,dZ_Rmor'],
                     ['Zoopl_E,dZ_Ea', 'Zoopl_E,dZ_Ec', 'Zoopl_E,dZ_Emor']]

    # color maps
    cmap = cmocean.cm.balance

    # center at zero?
    cmap_diverging = True

    # don't let there be a log scale for this one! it goes negative
    log_scale = False

elif rate_name=='oxycon':

    # title for figure
    rate_title = 'Oxygen Consumption'

    # grams of what in the units?   
    grams_of_what = 'O'

    # list of balance tables
    balance_table_list = ['oxy_Table.csv']

    # multiplier for each balance table
    multiplier_list = [-1,-1]

    # for each balance table, list of reactions to sum
    reaction_list = [['OXY,dOxCon','OXY,dMinDetCS1']]
    
    # color maps
    cmap = cmocean.cm.dense  

    # center at zero?
    cmap_diverging = False

    # include a log scale plot?
    log_scale = True

##################
### MAIN
##################

# check if plot directory exists and create it if not
if not os.path.exists(figure_path):
    os.makedirs(figure_path)

# time windows 
if time_period == 'Annual':
    time_titles = ['WY%d (Annual Average)' % wy]
    time_windows = [['%d-10-01' % (wy-1),'%d-10-01' % wy]]
elif time_period == 'Seasonal':
    time_titles = ['WY%d: Oct, Nov, Dec' % wy,
                   'WY%d: Jan, Feb, Mar' % wy,
                   'WY%d: Apr, May, Jun' % wy,
                   'WY%d: Jul, Aug, Sep' % wy]
    time_windows = [['%d-10-01' % (wy-1),'%d-01-01' % wy],
                    ['%d-01-01' % wy,'%d-04-01' % wy],
                    ['%d-04-01' % wy,'%d-07-01' % wy],
                    ['%d-07-01' % wy,'%d-10-01' % wy]]
elif time_period == 'Monthly':
    time_titles = ['WY%d: Oct' % wy,
                   'WY%d: Nov' % wy,
                   'WY%d: Dec' % wy,
                   'WY%d: Jan' % wy,
                   'WY%d: Feb' % wy,
                   'WY%d: Mar' % wy,
                   'WY%d: Apr' % wy,
                   'WY%d: May' % wy,
                   'WY%d: Jun' % wy,
                   'WY%d: Jul' % wy,
                   'WY%d: Aug' % wy,
                   'WY%d: Sep' % wy]
    time_windows = [['%d-10-01' % (wy-1),'%d-11-01' % (wy-1)],
                    ['%d-11-01' % (wy-1),'%d-12-01' % (wy-1)],
                    ['%d-12-01' % (wy-1),'%d-01-01' % wy],
                    ['%d-01-01' % wy,'%d-02-01' % wy],
                    ['%d-02-01' % wy,'%d-03-01' % wy],
                    ['%d-03-01' % wy,'%d-04-01' % wy],
                    ['%d-04-01' % wy,'%d-05-01' % wy],
                    ['%d-05-01' % wy,'%d-06-01' % wy],
                    ['%d-06-01' % wy,'%d-07-01' % wy],
                    ['%d-07-01' % wy,'%d-08-01' % wy],
                    ['%d-08-01' % wy,'%d-09-01' % wy],
                    ['%d-09-01' % wy,'%d-10-01' % wy]]
nwindows = len(time_titles)
assert nwindows==len(time_windows)

# units for normalized data
if norm == 'Volume':
    norm_units = 'g %s/m$^3$/d' % grams_of_what
elif norm == 'Area':
    norm_units = 'g %s/m$^2$/d' % grams_of_what
elif norm == 'None':
    norm_units = 'g %s/d' % grams_of_what
else:
    raise Exception("must specifiy normalization correctly")
norm_units_log = 'Log10 ' + norm_units

#####################
# here's the meat
#####################

# load the shapefile
gdf = gpd.read_file(shp_fn)
        
# get outline for plotting
gdf['Rate'] = 1
outline = gdf.iloc[iplot].dissolve(by='Rate')

# copy for plotting log
if log_scale:
    gdf_log = gdf.copy(deep=True)

# check number of balance tables 
ntables = len(balance_table_list)

# read the first balance table and sum up the reacitons, mutiplying by the appropriate stoichiometric multiplier
df = pd.read_csv(os.path.join(balance_table_path,balance_table_list[0]))
rate = multiplier_list[0]*df[reaction_list[0]].sum(axis=1)

# if there is more than one balance table, add those up too, summing the reactions by the multipliers
if ntables>1:
    for i in range(1,ntables):
        df = pd.read_csv(os.path.join(balance_table_path,balance_table_list[i]))
        rate = rate + multiplier_list[i]*df[reaction_list[i]].sum(axis=1)

# from the last balance table, grab all the control voulme ID and geometry info, as well as the time axis, then
# add the rate 
df = df[['time', 'Control Volume', 'Volume', 'Area']]
df['Rate'] = rate
df['time'] = df['time'].astype('datetime64[ns]')

# initiliaze color bar limits
pmin_all = 0
pmax_all = 0
if log_scale:
    pmin_log_all = 1e9
    pmax_log_all = -1e9

# initialize list of geodataframes
gdf_all = []
if log_scale:
    gdf_log_all = []

# loop through time windows, building up list of geodataframes, and finding the max and min values
for iw in range(nwindows):
    
    # get start and end time and boolean array for the time winodw
    time_window = time_windows[iw]
    time_start = np.datetime64(time_window[0])
    time_end = np.datetime64(time_window[1])
    timeall = df['time'].values
    indt = np.logical_and(timeall>=time_start,timeall<time_end)

    # loop through teh polygons
    for poly in range(len(gdf)):
    
        # polygon name
        polyname = 'polygon%d' % poly
        
        # get divisor based on normaliization
        if norm=='Volume':
            V = df.loc[df['Control Volume'] == polyname]['Volume'].iloc[0]
        elif norm=='Area':
            V = df.loc[df['Control Volume'] == polyname]['Area'].iloc[0]
        elif norm=='None':
            V = 1
            
        # find intersection of polygon and time window
        indp = (df['Control Volume'] == polyname).values
        ind = np.logical_and(indt,indp)
        
        # take the average over this time window, normalize, and load into a geodataframe
        gdf['Rate'].iloc[poly] = df.loc[ind]['Rate'].mean()/V
        if log_scale:
            gdf_log['Rate'].iloc[poly] = np.log10(df.loc[ind]['Rate'].mean()/V)
            

    # get the max and min parameter value
    pmax = np.nanpercentile(gdf['Rate'].iloc[iplot],cper)
    pmin = np.nanpercentile(gdf['Rate'].iloc[iplot],100-cper)
    if log_scale:
        pmax_log = np.nanpercentile(gdf_log['Rate'].iloc[iplot],cper)
        pmin_log = np.nanpercentile(gdf_log['Rate'].iloc[iplot],100-cper)
        
    
    # keep track of the biggest limits across time windows
    pmin_all = np.min([pmin,pmin_all])
    pmax_all = np.max([pmax,pmax_all])
    if log_scale:
        pmin_log_all = np.min([pmin_log,pmin_log_all])
        pmax_log_all = np.max([pmax_log,pmax_log_all])
    
    # append geodataframes
    gdf_all.append(copy.deepcopy(gdf))
    if log_scale:
        gdf_log_all.append(copy.deepcopy(gdf_log))

# if it's a diverging colormap, make the max and min the same amplitude
if cmap_diverging:
    pmax_all = np.max([pmax_all,-pmin_all])
    pmin_all = -pmax_all

# loop through time winodws again and plot
for iw in range(nwindows):
    
    # make figure title
    figure_title = '%s\n%s (%s)' % (time_titles[iw], rate_title, norm_units)
    if log_scale:
        figure_title_log = '%s\n%s (%s)' % (time_titles[iw], rate_title, norm_units_log)

    # make figure filename
    figure_fn = '%s_%s_%s_Map_%s_%04d.png' % (run_ID, rate_name, norm, time_period, iw)
    if log_scale:
        figure_log_fn = '%s_%s_log10_%s_Map_%s_%04d.png' % (run_ID, rate_name, norm, time_period, iw)
    
    # get geodataframe
    gdf = gdf_all[iw]
    if log_scale:
        gdf_log = gdf_log_all[iw]
    
    # open figure, append handles to master list, and turn axes off
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.axis('off')
    if log_scale:
        fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
        ax1.axis('off')

    # plot and save
    gdf.iloc[iplot].plot(ax=ax, column='Rate', cmap=cmap, vmin = pmin_all, vmax = pmax_all, legend=True)
    outline.boundary.plot(ax=ax,edgecolor='k')
    ax.set_title(figure_title)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_path, figure_fn))
    
    # log plot and save
    if log_scale:
        gdf_log.iloc[iplot].plot(ax=ax1, column='Rate', cmap=cmap, vmin = pmin_log_all, vmax = pmax_log_all, legend=True)
        outline.boundary.plot(ax=ax1,edgecolor='k')
        ax1.set_title(figure_title_log)
        fig1.tight_layout()
        fig1.savefig(os.path.join(figure_path, figure_log_fn))
    
    plt.close('all')