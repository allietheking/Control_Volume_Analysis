
''' allie'''

##################
# IMPORT MODULES
##################

import copy
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import os, sys
import pandas as pd
import cmocean 
import geopandas as gpd
    

###################
## USER INPUT
###################

# run ID list (can compare )
#run_ID_list = ['FR13_025','FR17_018','FR18_006']
run_ID_list = ['G141_13to18_207','G141_13to18_207','G141_13to18_207','G141_13to18_207','G141_13to18_207','G141_13to18_207']
#run_ID_list = ['G141_13to18_197','FR13_025','G141_13to18_197','G141_13to18_197','G141_13to18_197','G141_13to18_197','FR17_018','G141_13to18_197','FR18_006']

# water year
#wy = 2013
#wy = 2017
#wy = 2018
#wy_list = [2013, 2017, 2018]
wy_list = [2013, 2014, 2015, 2016, 2017, 2018]
#wy_list = [2013, 2013, 2014, 2015, 2016, 2017, 2017, 2018, 2018]

# specify the rates you want to plot right now (scroll down to see their definitions)
#rate_list = ['denit','dpp','dpp-benthic','dpp-pelagic','n-dpp','din_loss','din_recycling','dmin_water','dmin_sed','tn_loss','n-algae-sed','pon-sed','oxycon']
rate_list = ['dpp','dpp-benthic','dpp-pelagic']

# list of normalizations (divide by area, volume, or nothing)
#norm_list = ['None','Area','Volume']
norm_list = ['Area']

# list of time averaging periods (choices are annual, seasonal, monthly)
#time_period_list = ['Annual','Seasonal','Monthly']
time_period_list = ['Seasonal']

# figure size and subplots
figsize = (4*nruns,5)

# base directory (depends if running from laptop or on an hpc)
#base_dir = r'X:\hpcshared'
base_dir = '/chicagovol1/hpcshared'
base_dir_inputs = '/richmondvol1/hpcshared'

# path to balance tables (they should be organized by run ID within this folderr)
balance_table_path = os.path.join(base_dir,'NMS_Projects','Control_Volume_Analysis','Balance_Tables')

# nanpercentile for color map cutoff
cper = 97.5

# axis limits
axlim = (531280.627955355, 610527.2330959855, 4138850.8659710716, 4233404.6647333605)

# path to the shapefile for full res / aggregated runs
shp_fn_FR = os.path.join(base_dir_inputs,'inputs','shapefiles','Agg_mod_contiguous_plus_subembayments_shoal_channel.shp')
shp_fn_AGG = os.path.join(base_dir_inputs,'inputs','shapefiles','Agg_mod_contiguous_141.shp')

# subset of polygons to include in plot
iplot_FR = [0, 117, 139, 2, 113, 114, 115, 111, 1, 3, 116, 7, 4, 112, 5, 108, 109, 110, 9, 107, 
6, 8, 29, 10, 12, 11, 138, 137, 100, 19, 18, 13, 20, 26, 25, 21, 14, 24, 101, 102, 103, 104, 
105, 106, 28, 27, 22, 15, 30, 23, 35, 34, 33, 32, 31, 17, 41, 40, 39, 38, 37, 36, 16, 140, 
44, 43, 42, 46, 45, 49, 47, 97, 98, 86, 96, 95, 85, 93, 52, 87, 94, 48, 51, 50, 90, 88, 91, 
89, 92, 53, 56, 99, 54, 55, 57, 67, 65, 66, 136, 58, 59, 63, 64, 60, 62, 61, 143, 144, 141, 
68, 69, 70, 71, 78, 72, 84, 79, 146, 80, 82, 83, 142, 145, 73, 81, 74, 76, 77, 75]
iplot_AGG = [  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
    13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
    26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
    39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
    52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
    65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
    78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
    91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
   104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
   134, 136, 137, 138, 139, 140]

# output folder path
if all(np.array(run_ID_list)==run_ID_list[0]):
    figure_path = os.path.join(base_dir,'NMS_Projects','Control_Volume_Analysis','Plots',run_ID_list[0],'rate_maps')
else:
    figure_path = os.path.join(base_dir,'NMS_Projects','Control_Volume_Analysis','Plots','Compare_Water_Years','rate_maps')

##################
### MAIN
##################

# number of runs to compare
nruns = len(run_ID_list)
assert nruns == len(wy_list)

# check if plot directory exists and create it if not
if not os.path.exists(figure_path):
    os.makedirs(figure_path)

# what rate to plot? options are 
# 'denit','dpp','dpp-benthic','dpp-pelagic','n-dpp','din_loss','din_recycling','dmin_water','dmin_sed','tn_loss','n-algae-sed','pon-sed','oxycon'
for rate_name in rate_list:

    # time period (Annual, Seasonal, Monthly)
    for time_period in time_period_list:

        # how to norm
        for norm in norm_list:
    
            # here's where we set up what goes into these different rates
            if rate_name=='denit':
                
                # title for figure
                rate_title = 'Denitrification'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['din_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [["NO3,dDenit"]]
            
                # color maps
                cmap = cmocean.cm.dense
            
                # center at zero?
                cmap_diverging = False
            
                # include a log scale plot?
                log_scale = True
            
            elif rate_name=='dpp':
            
                # title for figure
                rate_title = 'Net Primary Productivity'
            
                # grams of what in the units?   
                grams_of_what = 'C'
            
                # list of balance tables
                balance_table_list = ['algae_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['Algae,dPPAlgae','Algae,dcPPAlgae']]
                
                # color maps
                cmap = cmocean.cm.algae
                #cmap = mpl.cm.seismic
            
                # center at zero?
                cmap_diverging = False
                #cmap_diverging = True
            
                # include a log scale plot?
                log_scale = True

            elif rate_name=='dpp-benthic':
            
                # title for figure
                rate_title = 'Net Primary Productivity (Benthic)'
            
                # grams of what in the units?   
                grams_of_what = 'C'
            
                # list of balance tables
                balance_table_list = ['diats1_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['DiatS1,dPPDiatS1']]
                
                # color maps
                cmap = cmocean.cm.algae
                #cmap = mpl.cm.seismic
            
                # center at zero?
                cmap_diverging = False
                #cmap_diverging = True
            
                # include a log scale plot?
                log_scale = True

            elif rate_name=='dpp-pelagic':
            
                # title for figure
                rate_title = 'Net Primary Productivity (Pelagic)'
            
                # grams of what in the units?   
                grams_of_what = 'C'
            
                # list of balance tables
                balance_table_list = ['diat_Table.csv','green_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1,1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['Diat,dPPDiat', 'Diat,dcPPDiat'],
                                 ['Green,dPPGreen', 'Green,dcPPGreen']]
                
                # color maps
                cmap = cmocean.cm.algae
                #cmap = mpl.cm.seismic
            
                # center at zero?
                cmap_diverging = False
                #cmap_diverging = True
            
                # include a log scale plot?
                log_scale = True
            
            elif rate_name=='n-dpp':
            
                # title for figure
                rate_title = 'Net Primary Productivity'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['n-algae_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['Algae,dPPAlgae','Algae,dcPPAlgae']]
                
                # color maps
                cmap = cmocean.cm.algae
            
                # center at zero?
                cmap_diverging = False
            
                # include a log scale plot?
                log_scale = True
            
            elif rate_name=='din_loss':
            
                # title for figure
                rate_title = 'DIN Assimilation'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['din_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum (this is all the reactions in the table)
                reaction_list = [['NH4,dMinDetN', 
                                  'DIN,dDINUpt', 
                                  'NO3,dDenit',
                                  'NO3,dNiDen', 
                                  'NH4,dMinPON1', 
                                  'NH4,dMinDON', 
                                  'NH4,dZ_NRes',
                                  'NH4,dNH4Aut',
                                  'ZERO: NH4,dNH4AUTS1', 
                                  'ZERO: NH4,dNitrif + NO3,dNitrif']]
            
                # color maps
                cmap = cmocean.cm.balance
            
                # center at zero?
                cmap_diverging = True
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = False

            elif rate_name=='din_recycling':
            
                # title for figure
                rate_title = 'DIN Recycling (Respiration + Mortality)'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['din_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum (this is all the reactions in the table)
                reaction_list = [["NH4,dZ_NRes",
                                  "NH4,dNH4Aut"]]
            
                # color maps
                cmap = cmocean.cm.balance
            
                # color maps
                cmap = cmocean.cm.amp
            
                # center at zero?
                cmap_diverging = False
            
                # include a log scale plot?
                log_scale = True

            elif rate_name=='dmin_water':
            
                # title for figure
                rate_title = 'DON + PON Mineralization'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['din_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum (this is all the reactions in the table)
                reaction_list = [["NH4,dMinPON1","NH4,dMinDON"]]
            
                # color maps
                cmap = cmocean.cm.matter
            
                # center at zero?
                cmap_diverging = False
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = True

            elif rate_name=='dmin_sed':
            
                # title for figure
                rate_title = 'Detritus Mineralization'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['din_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [1]
            
                # for each balance table, list of reactions to sum (this is all the reactions in the table)
                reaction_list = [["NH4,dMinDetN"]]
            
                # color maps
                cmap = cmocean.cm.turbid
            
                # center at zero?
                cmap_diverging = False
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = True

            elif rate_name=='tn_loss':
            
                # title for figure
                rate_title = 'TN Reactive Loss'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['tn_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['NO3,dDenit', 
                                  'NO3,dNiDen', 
                                  'Algae,dSedAlgae',
                                  'PON1,dSedPON1', 
                                  'DiatS1,dMrtDiatS1', 
                                  'DiatS1,dBurS1Diat',
                                  'NH4,dMinDetNS', 
                                  'NH4,dClam_NRes', 
                                  'PON1,dClam_NDef',
                                  'Algae,dClam_Algae', 
                                  'PON1,dClam_PON1' 
                                  #'ZERO: PON1,dCnvPPON1',
                                  #'ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN',
                                  #'ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet',
                                  #'ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor',
                                  #'ZERO: NH4,dNH4AUTS1 + DiatS1,dResS1Diat + DiatS1,dSWBuS1Dia + DiatS1,dDigS1Diat',
                                  #'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt',
                                  #'ZERO: NH4,dNH4UptS1 + NH4,dNH4US1D + NO3,dNO3UptS1 + DiatS1,dPPDiatS1',
                                  #'ZERO: NO3,dNitrif + NH4,dNitrif',
                                  #'ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1',
                                  #'ZERO: NH4,dMinPON1 + PON1,dMinPON1', 'ZERO: NH4,dMinDON + DON,dMinDON',
                                  #'ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut',
                                  #'ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec'
                                  ]]
                # color maps
                cmap = cmocean.cm.balance
            
                # center at zero?
                cmap_diverging = True
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = False



            elif rate_name=='n-algae-sed':
            
                # title for figure
                rate_title = 'Settling of Algae'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['tn_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [["Algae,dSedAlgae"
                                  #"NH4,dMinDetNS",
                                  #"NH4,dClam_NRes",
                                  #"PON1,dClam_NDef",
                                  #"Algae,dClam_Algae",
                                  #"PON1,dClam_PON1",
                                  #"ZERO: NO3,dNiDen",
                                  #"ZERO: PON1,dCnvPPON1",
                                  #"ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN",
                                  #"ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet",
                                  #"ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor",
                                  #"ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt",
                                  #"ZERO: NO3,dNitrif + NH4,dNitrif",
                                  #"ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1",
                                  #"ZERO: NH4,dMinPON1 + PON1,dMinPON1",
                                  #"ZERO: NH4,dMinDON + DON,dMinDON",
                                  #"ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut",
                                  #"ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec"]]
                                  ]]
                # color maps
                cmap = cmocean.cm.algae
            
                # center at zero?
                cmap_diverging = False
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = True

            elif rate_name=='pon-sed':
            
                # title for figure
                rate_title = 'Settling of PON'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['tn_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [[
                                  "PON1,dSedPON1"
                                  #"NH4,dMinDetNS",
                                  #"NH4,dClam_NRes",
                                  #"PON1,dClam_NDef",
                                  #"Algae,dClam_Algae",
                                  #"PON1,dClam_PON1",
                                  #"ZERO: NO3,dNiDen",
                                  #"ZERO: PON1,dCnvPPON1",
                                  #"ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN",
                                  #"ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet",
                                  #"ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor",
                                  #"ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt",
                                  #"ZERO: NO3,dNitrif + NH4,dNitrif",
                                  #"ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1",
                                  #"ZERO: NH4,dMinPON1 + PON1,dMinPON1",
                                  #"ZERO: NH4,dMinDON + DON,dMinDON",
                                  #"ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut",
                                  #"ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec"]]
                                  ]]
                # color maps
                cmap = mpl.cm.Greys
            
                # center at zero?
                cmap_diverging = False
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = True
            
            elif rate_name=='tn_include_sediment_loss':
            
                # title for figure
                rate_title = 'TN (including sediment N) Reactive Loss'
            
                # grams of what in the units?   
                grams_of_what = 'N'
            
                # list of balance tables
                balance_table_list = ['tn_include_sediment_Table.csv']
            
                # multiplier for each balance table
                multiplier_list = [-1]
            
                # for each balance table, list of reactions to sum
                reaction_list = [['NO3,dDenit', 
                                  'NO3,dNiDen',
                                  'DetNS2,dBurS2DetN', 
                                  'DiatS1,dBurS1Diat' 
                                  #'ZERO: PON1,dCnvPPON1',
                                  #'ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN',
                                  #'ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet',
                                  #'ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor',
                                  #'ZERO: DetNS1,dZ_NMrtS1 + DetNS1,dZ_DNS1 + DetNS1,dM_DNS1 + DetNS1,dG4_DNS1',
                                  #'ZERO: DetNS1,dSWMinDNS1 + DetNS1,dResS1DetN + DetNS1,dSWBuS1DtN + DetNS1,dDigS1DetN',
                                  #'ZERO: DetNS2,dSWMinDNS2 + DetNS2,dResS2DetN + DetNS2,dDigS1DetN + DetNS2,dDigS2DetN',
                                  #'ZERO: Mussel_R,dM_SpwDet + Grazer4_R,dG4_SpwDet',
                                  #'ZERO: NH4,dNH4AUTS1 + DiatS1,dResS1Diat + DiatS1,dSWBuS1Dia + DiatS1,dDigS1Diat',
                                  #'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt',
                                  #'ZERO: NH4,dNH4UptS1 + NH4,dNH4US1D + NO3,dNO3UptS1 + DiatS1,dPPDiatS1',
                                  #'ZERO: NO3,dNitrif + NH4,dNitrif',
                                  #'ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1',
                                  #'ZERO: NH4,dMinPON1 + PON1,dMinPON1', 'ZERO: NH4,dMinDON + DON,dMinDON',
                                  #'ZERO: DetNS1,dMrtDetNS1 + DiatS1,dMrtDiatS1',
                                  #'ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut',
                                  #'ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec',
                                  #'ZERO: PON1,dSedPON1 + DetNS1,dSedPON1',
                                  #'ZERO: NH4,dMinDetNS1 + NH4,dMinDetNS2 + DetNS1,dMinDetNS1 + DetNS2,dMinDetNS2',
                                  #'ZERO: Mussel_V,dM_Vmor + Mussel_E,dM_Emor + Mussel_R,dM_Rmor + DetNS1,dM_NMrtS1',
                                  #'ZERO: Grazer4_V,dG4_Vmor + Grazer4_E,dG4_Emor + Grazer4_R,dG4_Rmor + DetNS1,dG4_NMrtS1',
                                  #'ZERO: NH4,dM_NRes + PON1,dM_NDef + PON1,dM_PON1 + Diat,dM_Diat + Mussel_V,dM_Vgr + Mussel_E,dM_Ea + Mussel_E,dM_Ec + Mussel_R,dM_Rgr',
                                  #'ZERO: NH4,dG4_NRes + PON1,dG4_NDef + PON1,dG4_PON1 + Diat,dG4_Diat + Grazer4_V,dG4_Vgr + Grazer4_E,dG4_Ea + Grazer4_E,dG4_Ec + Grazer4_R,dG4_Rgr',
                                  #'ZERO: DetNS1,dBurS1DetN + DetNS2,dBurS1DetN',
                                  #'ZERO: Diat,dSedDiat + Green,dSedGreen + DetNS1,dSedAlgN'
                                  ]]
            
                # color maps
                cmap = mpl.cm.Spectral_r
            
                # center at zero?
                cmap_diverging = False
            
                # don't let there be a log scale for this one! it goes negative
                log_scale = True
            
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
            
            # time windows 
            wy_ref = 2009 # reference water year, pick something, doesn't matter what
            if time_period == 'Annual':
                time_titles = ['Annual Average']
                time_windows = [['%d-10-01' % (wy_ref-1),'%d-10-01' % wy_ref]]
            elif time_period == 'Seasonal':
                time_titles = ['Oct, Nov, Dec',
                               'Jan, Feb, Mar',
                               'Apr, May, Jun',
                               'Jul, Aug, Sep']
                time_windows = [['%d-10-01' % (wy_ref-1),'%d-01-01' % wy_ref],
                                ['%d-01-01' % wy_ref,'%d-04-01' % wy_ref],
                                ['%d-04-01' % wy_ref,'%d-07-01' % wy_ref],
                                ['%d-07-01' % wy_ref,'%d-10-01' % wy_ref]]
            elif time_period == 'Monthly':
                time_titles = ['Oct',
                               'Nov',
                               'Dec',
                               'Jan',
                               'Feb',
                               'Mar',
                               'Apr',
                               'May',
                               'Jun',
                               'Jul',
                               'Aug',
                               'Sep']
                time_windows = [['%d-10-01' % (wy_ref-1),'%d-11-01' % (wy_ref-1)],
                                ['%d-11-01' % (wy_ref-1),'%d-12-01' % (wy_ref-1)],
                                ['%d-12-01' % (wy_ref-1),'%d-01-01' % wy_ref],
                                ['%d-01-01' % wy_ref,'%d-02-01' % wy_ref],
                                ['%d-02-01' % wy_ref,'%d-03-01' % wy_ref],
                                ['%d-03-01' % wy_ref,'%d-04-01' % wy_ref],
                                ['%d-04-01' % wy_ref,'%d-05-01' % wy_ref],
                                ['%d-05-01' % wy_ref,'%d-06-01' % wy_ref],
                                ['%d-06-01' % wy_ref,'%d-07-01' % wy_ref],
                                ['%d-07-01' % wy_ref,'%d-08-01' % wy_ref],
                                ['%d-08-01' % wy_ref,'%d-09-01' % wy_ref],
                                ['%d-09-01' % wy_ref,'%d-10-01' % wy_ref]]
            
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
            
            # check number of balance tables 
            ntables = len(balance_table_list)
            
            # initiliaze color bar limits, which will be calculated by looking across all the runs and all the time windows
            pmin_all = 0
            pmax_all = 0
            if log_scale:
                pmin_log_all = 1e30
                pmax_log_all = -1e30
            
            # initialize dictionary containing geodataframes for plotting all the runs and time windows (this is confusing, sorry about that, but it works and is fast enough)
            gdf_all_dict = {}
            gdf_log_all_dict = {}
            
            # loop through the runs, building up a the dictionary of lists of geodataframes for plotting, and tracking the color bar limits 
            for irun in range(nruns):
            
                # get run ID and water year
                run_ID = run_ID_list[irun]
                wy = wy_list[irun]

                # path to the shapefile w/ the base level control volumes
                if 'FR' in run_ID:
                    shp_fn = shp_fn_FR
                    iplot = iplot_FR
                else:
                    shp_fn = shp_fn_AGG
                    iplot = iplot_AGG

                # load the shapefile and initialize Rate column
                gdf = gpd.read_file(shp_fn)
                gdf['Rate'] = np.nan
            
                # copy for plotting log
                if log_scale:
                    gdf_log = gdf.copy(deep=True)
            
                # read the first balance table and sum up the reacitons, mutiplying by the appropriate stoichiometric multiplier
                df = pd.read_csv(os.path.join(balance_table_path,run_ID,balance_table_list[0]))
                rate = multiplier_list[0]*df[reaction_list[0]].sum(axis=1)
                
                # if there is more than one balance table, add those up too, summing the reactions by the multipliers
                if ntables>1:
                    for i in range(1,ntables):
                        df = pd.read_csv(os.path.join(balance_table_path,run_ID,balance_table_list[i]))
                        rate = rate + multiplier_list[i]*df[reaction_list[i]].sum(axis=1)
                
                # from the last balance table, grab all the control voulme ID and geometry info, as well as the time axis, then
                # add the rate 
                df = df[['time', 'Control Volume', 'Volume', 'Area']]
                df['Rate'] = rate
                df['time'] = df['time'].astype('datetime64[ns]')
                
                # initialize list of geodataframes
                gdf_all = []
                if log_scale:
                    gdf_log_all = []
                
                # loop through time windows, building up list of geodataframes, and finding the max and min values
                for iw in range(nwindows):
                    
                    # get start and end time and boolean array for the time winodw
                    time_window = time_windows[iw]
                    time_start = np.datetime64(time_window[0]) - np.datetime64('%d-01-01' % wy_ref) + np.datetime64('%d-01-01' % wy)
                    time_end = np.datetime64(time_window[1]) - np.datetime64('%d-01-01' % wy_ref) + np.datetime64('%d-01-01' % wy)
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
                        
                    
                    # keep track of the biggest limits across time windows AND across all runs
                    pmin_all = np.min([pmin,pmin_all])
                    pmax_all = np.max([pmax,pmax_all])
                    if log_scale:
                        pmin_log_all = np.min([pmin_log,pmin_log_all])
                        pmax_log_all = np.max([pmax_log,pmax_log_all])
                    
                    # append geodataframes
                    gdf_all.append(copy.deepcopy(gdf))
                    if log_scale:
                        gdf_log_all.append(copy.deepcopy(gdf_log))
            
                # now append list of geodataframes to the dictionary
                gdf_all_dict[irun] = gdf_all.copy()
                if log_scale:
                    gdf_log_all_dict[irun] = gdf_log_all.copy()
                
            # if it's a diverging colormap, make the max and min the same amplitude
            if cmap_diverging:
                pmax_all = np.max([pmax_all,-pmin_all])
                pmin_all = -pmax_all
            
            # now loop through time winodws and plot all the runs we are comparing in the same figure with the same color bar limits
            for iw in range(nwindows):
                
                # make figure title
                cbar_title = '%s (%s)' % (rate_title, norm_units)
                if log_scale:
                    cbar_title_log = '%s (%s)' % (rate_title, norm_units_log)
            
                # make figure filename
                if all(run_ID_list[0]==np.array(run_ID_list)):
                    run_ID_string = run_ID + '_'
                else:
                    run_ID_string = ''
                    for run_ID in run_ID_list:
                        run_ID_string += run_ID + '_'
                figure_fn = '%s%s_%s_Map_%s_%04d.png' % (run_ID_string, rate_name, norm, time_period, iw)
                if log_scale:
                    figure_log_fn = '%s%s_log10_%s_Map_%s_%04d.png' % (run_ID_string, rate_name, norm, time_period, iw)
            
                # set up figure subwindows with room for a colorbar (looks best with 3 windows, hopefull ok if want other)
                fig = plt.figure(figsize=figsize)
                ax = ImageGrid(fig, 111,          # as in plt.subplot(111)
                             nrows_ncols=(1,nruns),
                             axes_pad=0.15,
                             share_all=True,
                             cbar_location="right",
                             cbar_mode="single",
                             cbar_size="7%",
                             cbar_pad=0.15,
                             ) 
                fig.suptitle(time_titles[iw])
            
                if log_scale: 
                    # open figure, append handles to master list, and turn axes off
                    fig1 = plt.figure(figsize=figsize)
                    ax1 = ImageGrid(fig1, 111,          # as in plt.subplot(111)
                                 nrows_ncols=(1,nruns),
                                 axes_pad=0.15,
                                 share_all=True,
                                 cbar_location="right",
                                 cbar_mode="single",
                                 cbar_size="7%",
                                 cbar_pad=0.15,
                                 ) 
                    fig1.suptitle(time_titles[iw])
            
                # loop through the runs and plot
                for irun in range(nruns):
            
                    # get geodataframe
                    gdf = gdf_all_dict[irun][iw]
                    if log_scale:
                        gdf_log = gdf_log_all_dict[irun][iw]

                    # path to the shapefile w/ the base level control volumes
                    if 'FR' in run_ID_list[irun]:
                        iplot = iplot_FR
                    else:
                        iplot = iplot_AGG

                    # get outline for plotting
                    gdf['dummy'] = 1
                    outline = gdf.iloc[iplot].dissolve(by='dummy')
                
                    # add to plot
                    gdf.iloc[iplot].plot(ax=ax[irun], column='Rate', cmap=cmap, vmin = pmin_all, vmax = pmax_all)
                    outline.boundary.plot(ax=ax[irun],edgecolor='k')
                    ax[irun].axis(axlim)
            
                    # turn off axis, set title 
                    ax[irun].axis('off')
                    ax[irun].set_title('WY%d (%s)' % (wy_list[irun],run_ID_list[irun]))
            
                    # same for log plot
                    if log_scale:
                        gdf_log.iloc[iplot].plot(ax=ax1[irun], column='Rate', cmap=cmap, vmin = pmin_log_all, vmax = pmax_log_all)
                        outline.boundary.plot(ax=ax1[irun],edgecolor='k')
                        ax1[irun].axis(axlim)
            
                        ax1[irun].axis('off')
                        ax1[irun].set_title('WY%d (%s)' % (wy_list[irun],run_ID_list[irun]))
            
                # add the colorbar
                norm1 = mpl.colors.Normalize(vmin=pmin_all, vmax=pmax_all)
                mpl.colorbar.ColorbarBase(ax[-1].cax, cmap=cmap,norm=norm1, label=cbar_title)
                if log_scale:
                    norm1 = mpl.colors.Normalize(vmin=pmin_log_all, vmax=pmax_log_all)
                    mpl.colorbar.ColorbarBase(ax1[-1].cax, cmap=cmap,norm=norm1, label=cbar_title_log)
            
                # save and close
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                fig.savefig(os.path.join(figure_path, figure_fn))
                if log_scale:
                    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
                    fig1.savefig(os.path.join(figure_path, figure_log_fn))
                
                plt.close('all')            