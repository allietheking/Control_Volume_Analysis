
'''
This script makes *.csv files with special variables that we output to the *.his file such as limnutdiat and rcgrodiat 
They are not really balance tables becuause we don't track the transport fluxes or mass, which for some of these variables
are meaningless anyhow because they are rates
'''

#################################################
# IMPORT MODULES (save stompy for later)
#################################################

import os, sys
import logging
import socket
import xarray as xr
import numpy as np
import pandas as pd
import datetime 
hostname = socket.gethostname()
if hostname == 'richmond':
    raise Exception("Do not run this script on richmond until we update the conda environment... run on chicago or your laptop instead")
try:
    import geopandas as gpd 
except:
    raise Exception('\nif geopandas is not found...\n' + 
                      'on chicago:  conda activate geo_env\n')


################
### USER INPUT
################

# run folder 
#run_folder = 'G141_13to18_146'
#run_folder = 'FR13_003' 
run_folder = 'FR13_025'
#run_folder = 'FR17_003'
#run_folder = 'FR17_018'
#run_folder = 'FR18_006'

# water year (for long agg run use WY13to18, this is part of path to the run files)
water_year = 'WY2013'
#water_year = 'WY2017'
#water_year = 'WY2018'
#water_year = 'WY13to18'

# delta run? (if so, script assumes it is full resolution)
is_delta = False

# is run full resolution? 
if 'FR' in run_folder:
    FR = True
else:
    FR = False

# root directory (for runs, shapefiles, and output csv files)
root_dir = '/richmondvol1/hpcshared/'   ## running on chicago
#root_dir = r'X:\hpcshared'              ## running on windows laptop with mounted drive

# stompy directory
stompy_dir = os.path.join(root_dir,'Software','stompy')

# output directory, including path
output_dir = os.path.join(root_dir, 'NMS_Projects','Control_Volume_Analysis','Balance_Tables')

# path to folder that contains run folder, which in turn contain model output including dwaq_hist.nc 
# and *-bal.his files. if there are no run subfolders, this is the direct path to dwaq_hist.nc and *-bal.his + 
# path to shapefiles defining polygons and transects for the monitoring regions (make sure these are the same ones used to
# initialize the DWAQ run)
if is_delta:
    shpfn_tran =  os.path.join(root_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Transects_Dave_Plus_WB_v4.shp') 
    shpfn_poly = os.path.join(root_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Polygons_Dave_Plus_WB_v4.shp')
    path = os.path.join(root_dir,'Delta','BGC_model','Full_res',run_folder)
    output_prefix = water_year.lower()
elif FR:
    shpfn_tran =  os.path.join(root_dir,'inputs','shapefiles','Agg_exchange_lines_plus_subembayments_shoal_channel.shp') # used in FR17003, FR13003
    shpfn_poly = os.path.join(root_dir,'inputs','shapefiles','Agg_mod_contiguous_plus_subembayments_shoal_channel.shp')
    path = os.path.join(root_dir,'Full_res','%s' % water_year)
    output_prefix = 'sfbay_dynamo000'
else:
    shpfn_poly =  os.path.join(root_dir,'inputs','shapefiles','Agg_mod_contiguous_141.shp')
    shpfn_tran = os.path.join(root_dir,'inputs','shapefiles','Agg_exchange_lines_141.shp')
    path = os.path.join(root_dir,'Grid141','%s' % water_year)
    output_prefix = 'sfbay_dynamo000'

# float format for csv files
float_format = '%1.16e'

# list of extra variables to put in the table
extra_variables =['limnutdiat',
 'limnutgree']
 #'limraddiat',
 #'limradgree',
 #'localdepth',
 #'rcgrodiat',
 #'rcgrogreen',
 #'rcrespdiat',
 #'rcrespgree',
 #'surf',
 #'totaldepth',
 #'depth',
 #'fppdiat',
 #'fppgreen',
 #'salinity',
 #'temp',
 #'volume']

#################
# IMPORT STOMPY
##################

if not stompy_dir in sys.path:
    sys.path.append(stompy_dir)
import stompy.model.delft.io as dio

##################
# MAIN
##################

# if delta, add "delta" prefix to the output folder
run_folder_out = run_folder
if is_delta:
    run_folder_out = 'Delta_' + run_folder_out

# if balance table output folder doesn't exist, create it
if not os.path.exists(os.path.join(output_dir,run_folder_out)):
    os.makedirs(os.path.join(output_dir,run_folder_out))
    
# path to his file
histfn      = os.path.join(path,run_folder,'dwaq_hist.nc')

# open hist file (may need to create it with stompy first)
hdata = xr.open_dataset(histfn)


# get the start time and make sure his and his-bal start times match
start_time = pd.to_datetime(hdata.time.values[0])
start_date = np.datetime64('%d-%02d-%02d' % (start_time.year, start_time.month, start_time.day))

# if the simulation does not start at midnight, subtract the time of day, and add it back later
offset_time = start_time.to_datetime64() - start_date
hdata['time'] = hdata['time'] - offset_time

# renumber the polygons so they match the shape file -- 
# newer version of stompy scrambles the numbers but does not scramble the order
polyc = 0
for i in range(len(hdata.location_names.values[0])):
    if 'polygon' in hdata.location_names.values[0][i]:
        hdata.location_names.values[0][i]='polygon%d' % polyc
        polyc = polyc + 1 

# create the balance tables directory if it doesn't exist yet
if not os.path.exists(os.path.join(output_dir,run_folder_out)):
    os.makedirs(os.path.join(output_dir,run_folder_out))

# create name for balance table output file
outfile = os.path.join(output_dir,run_folder_out,'Extra_variables_Table.csv')

# loop through all the parameters 
is_first = True
for varname in extra_variables:

    print(varname)
    
    ##%% Get all the data
    PolygonBL = ['polygon' in name for name in hdata.location_names.values[0]]
    indP = np.where(PolygonBL)[0]

    sys.exit()
    varP = hdata.isel(nSegment=indP)[varname]
    polygon_list = hdata.location_names.values[0][PolygonBL]

    # check the frequency of the data, and if frequency is higher than daily, resample onto a daily axis
    deltat_P = (varP.time[1]-varP.time[0]).values

    # ... if polygon output is less than daily, resample onto daily axis (take instantaneous snapshots)
    if deltat_P < np.timedelta64(1,'D'):
        varP = varP.resample(time='1D').nearest()
    elif deltat_P==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his file has time step greater than one day')
    
    # load shapefiles
    poly_df = gpd.read_file(shpfn_poly)
    pmax = len(poly_df)
    
    #%% collect variable for each polygon in a dataframe 
    df_var = pd.DataFrame()  
    for i,p in enumerate(polygon_list):   

        df = pd.DataFrame()
        df[varname] = varP.isel(nSegment=i).to_pandas()   
        df.insert(0, 'Control Volume', p)

        df_var = pd.concat([df_var, df]) 

    # add to dataframe with all the variable
    if is_first:
        df_output = df_var.copy()
        is_first = False
    else:
        df_output[varname] = df_var[varname].copy()    

sys.exit()


    
# adjust for offset time
time = pd.to_datetime(df_output.index + offset_time)
df_output.set_index(time, inplace=True)

# save
df_output.to_csv(outfile,columns=column_list,float_format=float_format)   
