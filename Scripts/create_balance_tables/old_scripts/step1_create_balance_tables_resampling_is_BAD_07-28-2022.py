
'''
This script converts *.his and *-bal.his data to *.csv formatted balance tables containing
daily fluxes and reaction terms in the "monitoring regions" defined by polygons and transects.
Updated by Allie in 2022 to run in new python environment on chicago:
    source activate geo_env
    cd /richmondvol1/hpcshared/NMS_Projects/Control_Volume/Scripts/create_balance_tables
and run from there
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
run_folder = 'G141_13to18_182'
#run_folder = 'FR13_003' 
#run_folder = 'FR13_025'
#run_folder = 'FR17_003'
#run_folder = 'FR17_018'
#run_folder = 'FR18_006'

# water year (for long agg run use WY13to18, this is part of path to the run files)
#water_year = 'WY2013'
#water_year = 'WY2017'
#water_year = 'WY2018'
water_year = 'WY13to18'

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

#################
# IMPORT STOMPY
##################

if not stompy_dir in sys.path:
    sys.path.append(stompy_dir)
import stompy.model.delft.io as dio

##################
# FUNCTIONS
##################

# define function
def Poly2Transect(left,right,pi='all'):
    # find the transects for each polygon and the signs based on the following
    # rule: transects with their 'from' segment left of the path are positive
    # otherwise negated; so left polygons are given a negative sign
    def p2t_i(i):
        indl = np.where(left==i)[0]
        indr = np.where(right==i)[0]
        signl = np.ones_like(indl)*-1
        signr = np.ones_like(indr)
        adj_poly_l = right[indl].values
        adj_poly_r = left[indr].values
        return {'transect':np.concatenate([indl,indr]),
                 'sign':np.concatenate([signl,signr]),
                 'adjacent': np.concatenate([adj_poly_l,adj_poly_r])}
    if pi=='all':
        p2t = []
        for i in np.arange(pmax):
            p2t.append(p2t_i(i))
    elif isinstance(pi,int):
        p2t = p2t_i(pi)
    else:
        raise ValueError("The type of pi is not implemented")             
    return p2t

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

# setup logging to file and to screen 
try:
    logging.basicConfig(
    level=logging.INFO,
    mode='w',
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(output_dir,run_folder_out,"log_step1.log")),
        logging.StreamHandler(sys.stdout)
    ])
except:
    logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(output_dir,run_folder_out,"log_step1.log")),
        logging.StreamHandler(sys.stdout)
    ])

# add some basic info to log file
user = os.getlogin()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('These balance tables were produced on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))
    
# path to his and his bal files
histfn      = os.path.join(path,run_folder,'dwaq_hist.nc')
histbal_fn  = os.path.join(path,run_folder,'dwaq_hist_bal.nc')

# log start of readin files
logging.info('reading %s and %s' % (histfn,histbal_fn))

# open hist file (may need to create it with stompy first)
hdata = xr.open_dataset(histfn)

# open hist-bal file, create if needed
try:
    hbdata = xr.open_dataset(histbal_fn)
except:
    fn = os.path.join(path,run_folder,'%s-bal.his' % output_prefix)
    hbdata = dio.bal_his_file_xarray(fn)
    hbdata.to_netcdf(histbal_fn)
    hbdata = xr.open_dataset(histbal_fn)

# get the start time and make sure his and his-bal start times match
start_time = pd.to_datetime(hdata.time.values[0])
if not start_time == pd.to_datetime(hbdata.time.values[0]):
    raise Exception('start time of his-bal data doesn\'t match start time of his data')
start_date = np.datetime64('%d-%02d-%02d' % (start_time.year, start_time.month, start_time.day))

# if the simulation does not start at midnight, subtract the time of day, and add it back later
offset_time = start_time.to_datetime64() - start_date
hdata['time'] = hdata['time'] - offset_time
hbdata['time'] = hbdata['time'] - offset_time

# renumber the polygons and transects so they match the shape file -- 
# newer version of stompy scrambles the numbers but does not scramble the order
polyc = 0
tranc = 0
for i in range(len(hdata.location_names.values[0])):
    if 'transect' in hdata.location_names.values[0][i]:
        hdata.location_names.values[0][i]='transect%04d' % tranc
        tranc = tranc + 1
    elif 'polygon' in hdata.location_names.values[0][i]:
        hdata.location_names.values[0][i]='polygon%d' % polyc
        polyc = polyc + 1 
polyc = 0
for i in range(len(hbdata.region.values)):
    if 'polygon' in hbdata.region.values[i]:
        hbdata.region.values[i] = 'polygon%d' % polyc
        polyc = polyc + 1

# create the balance tables directory if it doesn't exist yet
if not os.path.exists(os.path.join(output_dir,run_folder_out)):
    os.makedirs(os.path.join(output_dir,run_folder_out))
    
# loop through all the parameters (nh4, no3, diat, etc.)
varnames = [var.lower() for var in hbdata.sub.values]
for varname in varnames:
    
    # create name for balance table output file
    outfile = os.path.join(output_dir,run_folder_out,varname +'_Table.csv')
    
    ##%% Get all the data
    TransectBL = ['transect' in name for name in hdata.location_names.values[0]]
    indT = np.where(TransectBL)[0]
    PolygonBL = ['polygon' in name for name in hdata.location_names.values[0]]
    indP = np.where(PolygonBL)[0]
    PolygonBL_bal = ['polygon' in name for name in hbdata.region.values]
    indP_bal = np.where(PolygonBL_bal)[0]
    varT = hdata.isel(nSegment=indT)[varname]
    varP = hdata.isel(nSegment=indP)[varname]
    fieldBL = [varname.lower()+',' in name.lower() for name in hbdata.field.values]
    indF = np.where(fieldBL)[0]
    varP_bal = hbdata.isel(region=indP_bal).isel(field=indF)
    Vp = hdata.isel(nSegment=indP)['volume']  

    # check the frequency of the data, and if frequency is higher than daily, resample onto a daily axis
    deltat_P = (varP.time[1]-varP.time[0]).values
    deltat_T = (varT.time[1]-varT.time[0]).values
    deltat_B = (varP_bal.time[1]-varP_bal.time[0]).values

    # ... if polygon output is less than daily, resample onto daily axis (take instantaneous snapshots)
    if deltat_P < np.timedelta64(1,'D'):
        varP = varP.resample(time='1D').nearest()
        Vp_mean = Vp.resample(time='1D').mean()
        Vp = Vp.resample(time='1D').nearest()
    elif deltat_P==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his file has time step greater than one day')
    # ... transect output should be integrated in time, bizarrely the integral is open on the 
    # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
    # ... so add a day to it
    if deltat_T<np.timedelta64(1,'D'):
        varT = varT.resample(time='1D',closed='right').sum(dim='time')[0:-1,:]
        varT['time'] = varT['time'] + np.timedelta64(1,'D')
    elif deltat_T==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his file has time step greater than one day')
    # ... balance output should be integrated in time, bizarrely the integral is open on the 
    # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
    # ... so add a day to it
    if deltat_B<np.timedelta64(1,'D'):
        #varP_bal = varP_bal.resample(time='1D').sum(dim='time')[0:-1,:] # script stopped working on 6/29/2022, this fixed it
        varP_bal = varP_bal.resample(time='1D').sum(dim='time') # script stopped working on 6/29/2022, this fixed it
        #varP_bal['time'] = varP_bal['time'] + np.timedelta64(1,'D') # script stopped working on 6/29/2022, this fixed it
    elif deltat_B==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his-bal file has time step greater than one day')

    # now make sure polygon, transect, and balance data have same number of time steps (this 
    # condition may be violated in case of incomplete simulation)
    tmin = np.min([varP.time.values[-1],varT.time.values[-1],varP_bal.time.values[-1]])
    varP = varP.where(varP.time<=tmin,drop=True)
    varT = varT.where(varT.time<=tmin,drop=True)
    varP_bal = varP_bal.where(varP_bal.time<=tmin,drop=True)
    
    # load shapefiles
    gdf     = gpd.read_file(shpfn_tran)
    left    = gdf.left
    right   = gdf.right
    #pmax = max(left.max(),right.max()) +1 # total number of polygons and they are zero-based
    
    poly_df = gpd.read_file(shpfn_poly)
    pmax = len(poly_df)

    Area = poly_df.area.values
    
    #%% Getting the transect fluxes
    p2t = Poly2Transect(left,right)    
    
    #%% Outputting variables for each polygon. 
    diffVar = (varP*Vp).diff(dim='time')  
    varTv = varT.values                         
    df_output = pd.DataFrame()  
    for i,p in enumerate(varP_bal.region.values):
    
        pi_df = varP_bal.bal.sel(region=p).to_pandas() 
        pi_df['dVar/dt'] = np.append(diffVar[0,i],diffVar[:,i])
        
        if i==0:
            column_list = ['Control Volume','Concentration (mg/l)','Volume','Volume (Mean)',
                           'Area'] + list(pi_df.columns)     
        
        p2t_i = p2t[i]
        Fluxes = varTv[:,p2t_i['transect']]*p2t_i['sign']
        for t in np.arange(np.shape(Fluxes)[1]):
            cname = 'To_poly'+str(t)
            fname = 'Flux'+str(t)
            pi_df[cname] = p2t_i['adjacent'][t]
            pi_df[fname] = Fluxes[:,t]
            if cname not in column_list:
                column_list += [cname,fname]
        To_transect = np.nonzero(['To_poly' in name for name in pi_df.columns])[0]
        Others = np.nonzero(['To_poly' not in name for name in pi_df.columns])[0]
        pi_df_sum = pi_df.iloc[:,Others]
        pi_df_daily = pi_df.iloc[:,To_transect]
        pi_df_comb = pd.concat([pi_df_sum,pi_df_daily],axis=1)
        
        pi_df_comb['Concentration (mg/l)'] = varP.isel(nSegment=i).values
        pi_df_comb['Control Volume'] = p
        pi_df_comb['Volume'] = Vp.values[:,i]
        pi_df_comb['Volume (Mean)'] = Vp_mean.values[:,i] 
        pi_df_comb['Area'] = Area[i]    

        df_output = pd.concat([df_output,pi_df_comb])
            
    df_output = df_output[column_list]    

    
    # adjust for offset time
    time = pd.to_datetime(df_output.index + offset_time)
    df_output.set_index(time, inplace=True)

    # save
    df_output.to_csv(outfile,columns=column_list,float_format=float_format)   
    logging.info('Saved %s' % outfile)
