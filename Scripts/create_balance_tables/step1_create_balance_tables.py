
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
import xarray as xr
import numpy as np
import pandas as pd
import datetime 
import socket
hostname = socket.gethostname()
if hostname == 'richmond':
    raise Exception("Do not run this script on richmond until we update the conda environment... run on chicago or your laptop instead")
try:
    import geopandas as gpd 
except:
    raise Exception('\nif geopandas is not found...\n' + 
                      'on chicago:  conda activate geo_env\n')

# if running the script alone, load the configuration module (in this folder)
if __name__ == "__main__":

    import importlib
    import step0_config
    importlib.reload(step0_config)

################
### USER INPUT
################

# get variables out of the configuration module (see step0_config.py in this folder)
from step0_config import runid, is_delta, substance_list, balance_table_dir, model_inout_dir, stompy_dir, float_format

# for some parameters, the dwaq_hist.nc file does not have the correct units ... use this to override the units
units_override = {'zoopl_e' : 'gC/m3',
                  'zoopl_r' : 'gC/m3',
                  'zoopl_v' : 'gC/m3',
                  'zoopl_n' : '#/m3',}

# is run full resolution? 
if 'FR' in runid:
    FR = True
else:
    FR = False

# ugly complicated process to get the water year from the runid -- this is just for identifying the folder the run is stored in
if 'FR' in runid:
    # this extracts the 2 digit water year, assuming format of runid is like FR13_003 for WY2013 run 003
    yr = int(runid.split('_')[0][2:])
    # turn into water year string
    water_year = 'WY%d' % (2000 + yr)
elif 'G141' in runid:
    # there are two formats for agg runs, G141_13_003 is water year 2013, G141_13to18_207 is water years 2013-2018
    # get the string that represents the water year
    yr = runid.split('_')[1]
    # if it spans mutlple water years, keep the string, just add 'WY' in front of it
    if 'to' in yr:
        water_year = 'WY' + yr
    # otherwise extract the integer and add it to 2000
    else:
        water_year = 'WY%d' % (2000 + int(yr))

# path to folder that contains run folder, which in turn contain model output including dwaq_hist.nc 
# and *-bal.his files. if there are no run subfolders, this is the direct path to dwaq_hist.nc and *-bal.his + 
# path to shapefiles defining polygons and transects for the monitoring regions (make sure these are the same ones used to
# initialize the DWAQ run)
if is_delta:
    shpfn_tran =  os.path.join(model_inout_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Transects_Dave_Plus_WB_v4.shp') 
    shpfn_poly = os.path.join(model_inout_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Polygons_Dave_Plus_WB_v4.shp')
    path = os.path.join(model_inout_dir,'Delta','BGC_model','Full_res',runid)
    output_prefix = water_year.lower()
elif FR:
    shpfn_tran =  os.path.join(model_inout_dir,'inputs','shapefiles','Agg_exchange_lines_plus_subembayments_shoal_channel.shp') # used in FR17003, FR13003
    shpfn_poly = os.path.join(model_inout_dir,'inputs','shapefiles','Agg_mod_contiguous_plus_subembayments_shoal_channel.shp')
    path = os.path.join(model_inout_dir,'Full_res','%s' % water_year)
    output_prefix = 'sfbay_dynamo000'
else:
    shpfn_poly =  os.path.join(model_inout_dir,'inputs','shapefiles','Agg_mod_contiguous_141.shp')
    shpfn_tran = os.path.join(model_inout_dir,'inputs','shapefiles','Agg_exchange_lines_141.shp')
    path = os.path.join(model_inout_dir,'Grid141','%s' % water_year)
    output_prefix = 'sfbay_dynamo000'

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

def logger_cleanup():

    ''' check for open log files, close them, release handlers'''

    # clean up logging
    logger = logging.getLogger()
    handlers = list(logger.handlers)
    if len(handlers)>0:
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler) 


##################
# MAIN
##################

# if balance table output folder doesn't exist, create it
if not os.path.exists(balance_table_dir):
    os.makedirs(balance_table_dir)

# setup logging to file and to screen 
logger_cleanup()
logging.basicConfig(
level=logging.INFO,
format="%(asctime)s [%(levelname)s] %(message)s",
handlers=[
    logging.FileHandler(os.path.join(balance_table_dir,"log_step1.log"),'w'),
    logging.StreamHandler(sys.stdout)
])

# add some basic info to log file
user = os.getlogin()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('These balance tables were produced on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))
    
# log configuration variables
logging.info('The following global variables were loaded from step0_conf.py:')
logging.info('    runid = %s' % runid)
logging.info('    is_delta = %r' % is_delta)
logging.info('    substance_list = %s' % substance_list)
logging.info('    balance_table_dir = %s' % balance_table_dir)
logging.info('    model_inout_dir = %s' % model_inout_dir)
logging.info('    stompy_dir = %s' % stompy_dir)
logging.info('    float_format = %s' % float_format)

# path to his and his bal files
histfn      = os.path.join(path,runid,'dwaq_hist.nc')
histbal_fn  = os.path.join(path,runid,'dwaq_hist_bal.nc')

# log start of readin files
logging.info('reading %s and %s' % (histfn,histbal_fn))

# open hist file (may need to create it with stompy first)
hdata = xr.open_dataset(histfn)

# open hist-bal file, create if needed
try:
    hbdata = xr.open_dataset(histbal_fn)
except:
    fn = os.path.join(path,runid,'%s-bal.his' % output_prefix)
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
if not os.path.exists(balance_table_dir):
    os.makedirs(balance_table_dir)
    
# loop through all the parameters (nh4, no3, diat, etc.)
varnames = [var.lower() for var in hbdata.sub.values]
for varname in varnames:

    # if not processing all substances skip substances that aren't in the list
    if not substance_list=='all':
        if not varname in substance_list:
            logging.info('Skipping substance %s because it is not in the user-specified list of substances to process ...' % varname)     
            continue

    # otherwise create balance table
    logging.info('Creating balance table for substance %s...' % varname)
    
    # create name for balance table output file
    outfile = os.path.join(balance_table_dir,varname +'_Table.csv')
    
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
    # ... for volume additionally provide a mean value for normalizing rates
    if deltat_P < np.timedelta64(1,'D'):
        varP = varP.resample(time='1D').nearest()
        Vp = Vp.resample(time='1D').nearest()
        Vp_mean = Vp.resample(time='1D',closed='right',label='right').mean(dim='time')
    elif deltat_P==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his file has time step greater than one day')
    # ... transect output should be integrated in time, bizarrely the integral is open on the 
    # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
    # ... so add a day to it
    if deltat_T<np.timedelta64(1,'D'):
        varT = varT.resample(time='1D',closed='right',label='right').sum(dim='time')
    elif deltat_T==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his file has time step greater than one day')
    # ... balance output should be integrated in time, bizarrely the integral is open on the 
    # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
    # ... so add a day to it
    if deltat_B<np.timedelta64(1,'D'):
        varP_bal = varP_bal.resample(time='1D',closed='right',label='right').sum(dim='time')
    elif deltat_B==np.timedelta64(1,'D'):
        pass
    else:
        raise Exception('ERROR: his-bal file has time step greater than one day')

    # subtract one day from all output, becasue currently the values represent the backwards average 
    # over the PREVIOUS day
    varP['time'] = varP['time'] - np.timedelta64(1,'D')
    Vp['time'] = Vp['time'] - np.timedelta64(1,'D')
    Vp_mean['time'] = Vp_mean['time'] - np.timedelta64(1,'D')
    varT['time'] = varT['time'] - np.timedelta64(1,'D')
    varP_bal['time'] = varP_bal['time'] - np.timedelta64(1,'D')

    # now make sure polygon, transect, and balance data have same number of time steps (this 
    # condition may be violated in case of incomplete simulation)
    tmin = np.min([varP.time.values[-1],varT.time.values[-1],varP_bal.time.values[-1]])
    varP = varP.where(varP.time<=tmin,drop=True)
    Vp = Vp.where(Vp.time<=tmin,drop=True)
    Vp_mean = Vp_mean.where(Vp_mean.time<=tmin,drop=True)
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
    
    # get the units and determine if per area or per volume, then compute dMass/dt accordingly
    units1 = varP.units
    if varname in units_override.keys():
        units = units_override[varname]
        logging.info('units in dwaq_hist.nc are %s, overriding with %s' % (units1, units))
    else:
        units = units1
    if '/m2' in units:
        logging.info('units are %s, multiplying by area to get dVar/dt and converting concentration to volumetric' % units)
        diffVar = (varP*Area).diff(dim='time')
        Conc = (varP*Area)/Vp
    elif '/m3' in units:
        logging.info('units are %s, multiplying by volume to get dVar/dt' % units)
        diffVar = (varP*Vp).diff(dim='time') 
        Conc = varP.copy(deep=True)

    #%% Outputting variables for each polygon.  
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
        
        pi_df_comb['Concentration (mg/l)'] = Conc.isel(nSegment=i).values
        pi_df_comb['Control Volume'] = p
        pi_df_comb['Volume'] = Vp.values[:,i]
        pi_df_comb['Volume (Mean)'] = Vp_mean.values[:,i] 
        pi_df_comb['Area'] = Area[i]    

        # remove the first time step because it is garbage
        pi_df_comb = pi_df_comb.iloc[1:]

        df_output = pd.concat([df_output,pi_df_comb])
            
    df_output = df_output[column_list]    

    
    # adjust for offset time
    time = pd.to_datetime(df_output.index + offset_time)
    df_output.set_index(time, inplace=True)

    # save
    df_output.to_csv(outfile,columns=column_list,float_format=float_format)   
    logging.info('Saved %s' % outfile)

# clean up logging
logger_cleanup()