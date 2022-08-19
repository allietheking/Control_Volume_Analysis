
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
import matplotlib.pylab as plt
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

# monitoring stations
station_csv = os.path.join(root_dir,'inputs','stations','stations_v4.csv')

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

# read shapefile with polygons
poly_df = gpd.read_file(shpfn_poly)
pmax = len(poly_df)
poly_df['centroid'] = poly_df.centroid
xc_poly = []
yc_poly = []
for point in poly_df.centroid:
    xy = point.xy
    xc_poly.append(xy[0][0])
    yc_poly.append(xy[1][0])
poly_df['xc'] = xc_poly
poly_df['yc'] = yc_poly

# read monitoring station csv
station_df = pd.read_csv(station_csv).iloc[0:40]
smax = len(station_df)

# plot
fig,ax = plt.subplots(1,2,figsize=(14,7))
axlim = (577158.2467631227, 585504.1578651706, 4144651.964346677, 4152830.585471178)
ax1 = ax[0]
#axlim =  (484008.604365, 621165.308335, 4134823.9835110214, 4227288.720594713) # whole bay
poly_df.iloc[0:141].plot(ax=ax1, edgecolor='k')
ax1.axis(axlim)
for i in range(pmax):
    xc1 = poly_df.iloc[i].xc
    yc1 = poly_df.iloc[i].yc
    if xc1>=axlim[0] and xc1<=axlim[1] and yc1>=axlim[2] and yc1<=axlim[3]:
        ax1.text(xc1, yc1, 'P%d' % i)
for i in range(smax):
    xc1 = station_df.iloc[i]['utm_x']
    yc1 = station_df.iloc[i]['utm_y']
    if xc1>=axlim[0] and xc1<=axlim[1] and yc1>=axlim[2] and yc1<=axlim[3]:
        ax1.plot(xc1, yc1,'bo')
        ax1.text(xc1, yc1, station_df.iloc[i]['Station'])
ax1.axis('off')
axlim = (556364.286203918, 575903.4985746137, 4156455.8762679338, 4171147.1726718647)
ax1.set_title('Compare USGS Staiton 35 with Polygon 1')
ax1 = ax[1]
#axlim =  (484008.604365, 621165.308335, 4134823.9835110214, 4227288.720594713) # whole bay
poly_df.iloc[0:141].plot(ax=ax1, edgecolor='k')
ax1.axis(axlim)
for i in range(pmax):
    xc1 = poly_df.iloc[i].xc
    yc1 = poly_df.iloc[i].yc
    if xc1>=axlim[0] and xc1<=axlim[1] and yc1>=axlim[2] and yc1<=axlim[3]:
        ax1.text(xc1, yc1, 'P%d' % i)
for i in range(smax):
    xc1 = station_df.iloc[i]['utm_x']
    yc1 = station_df.iloc[i]['utm_y']
    if xc1>=axlim[0] and xc1<=axlim[1] and yc1>=axlim[2] and yc1<=axlim[3]:
        ax1.plot(xc1, yc1,'bo')
        ax1.text(xc1, yc1-500, station_df.iloc[i]['Station'])
ax1.axis('off')
ax1.set_title('Compare USGS Staiton 28 with Polygon 104')
plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Map_Station_Polygon_Comparison.png')
plt.close()

# compare station 28 and polygon 104 because their station 28 is at the centroid of polygon 104, and it's the channel, so hopefully kind of uniform???
# also compare station 35 and polygon 1 (in lower south bay) even though they don't overlap perfectly, to get a different depth

# time
time_bal = hbdata.time.values
time = hdata.time.values

# need to convert units of his-bal rates from mass per time step to mass per day
tmult_bal = np.timedelta64(1,'D')/(time_bal[1]-time_bal[0]) 


##%% Get all the balance dpp values for polygons and the stations
p104_bal = {}
PolygonBL_bal = 'polygon104' ==hbdata.region.values
indP_bal = np.where(PolygonBL_bal)[0]
for var in ['Diat,dPPDiat','Diat,dcPPDiat']:
    fieldBL = hbdata.field.values == var
    indF = np.where(fieldBL)[0]
    p104_bal[var] = hbdata.isel(region=indP_bal).isel(field=indF).bal[:,0,0].values * tmult_bal
p001_bal = {}
PolygonBL_bal = 'polygon1'==hbdata.region.values
indP_bal = np.where(PolygonBL_bal)[0]
for var in ['Diat,dPPDiat','Diat,dcPPDiat']:
    fieldBL = hbdata.field.values == var
    indF = np.where(fieldBL)[0]
    p001_bal[var] = hbdata.isel(region=indP_bal).isel(field=indF).bal[:,0,0].values * tmult_bal
s028_bal = {}
PolygonBL_bal = '28'==hbdata.region.values
indP_bal = np.where(PolygonBL_bal)[0]
for var in ['Diat,dPPDiat','Diat,dcPPDiat']:
    fieldBL = hbdata.field.values == var
    indF = np.where(fieldBL)[0]
    s028_bal[var] = hbdata.isel(region=indP_bal).isel(field=indF).bal[:,0,0].values * tmult_bal
s035_bal = {}
PolygonBL_bal = '35' ==hbdata.region.values
indP_bal = np.where(PolygonBL_bal)[0]
for var in ['Diat,dPPDiat','Diat,dcPPDiat']:
    fieldBL = hbdata.field.values == var
    indF = np.where(fieldBL)[0]
    s035_bal[var] = hbdata.isel(region=indP_bal).isel(field=indF).bal[:,0,0].values * tmult_bal


# list of variables
var_list = ['diat','volume','surf','totaldepth','depth','localdepth','eta','limnutdiat','limraddiat',
            'fppdiat','rcgrodiat','rcrespdiat','tavg_dppdiat','tavg_rcgrodiat','tavg_rcrespdiat']


# polygon data from hist
PolygonBL ='polygon104'==hdata.location_names.values[0]
indP = np.where(PolygonBL)[0][0]
p104 = {}
for var in var_list:
    p104[var] = hdata[var][:,indP]
    p104[var] = hdata[var][:,indP]
PolygonBL = 'polygon1'==hdata.location_names.values[0]
indP = np.where(PolygonBL)[0][0]
p001 = {}
for var in var_list:
    p001[var] = hdata[var][:,indP]

## overlapping usgs crusise station data
StationBL = '28                  '==hdata.location_names.values[0]
indP = np.where(StationBL)[0]
s028 = {}
for var in var_list:
    s028[var] = hdata[var][:,indP]
StationBL = '35                  '== hdata.location_names.values[0]
indP = np.where(StationBL)[0]
s035 = {}
for var in var_list:
    s035[var] = hdata[var][:,indP]

## overlapping usgs crusise station data, by layer
nlay = 10
ntim = len(time)
s028_lay = {}
for var in var_list:
    s028_lay[var] = np.zeros((ntim,nlay))
    for ilay in range(nlay):
        StationBL1 = ('28_%d                ' % ilay)== hdata.location_names.values[0]
        indP = np.where(StationBL1)[0]
        s028_lay[var][:,ilay] = hdata[var].values[:,indP][:,0]
s035_lay = {}
for var in var_list:
    s035_lay[var] = np.zeros((ntim,nlay))
    for ilay in range(nlay):
        StationBL1 = ('35_%d                ' % ilay)== hdata.location_names.values[0]
        indP = np.where(StationBL1)[0]
        s035_lay[var][:,ilay] = hdata[var].values[:,indP][:,0]

# get time axes for computing tavg
time_bal_tavg = np.unique(time_bal.astype('datetime64[D]'))
time_tavg = np.unique(time.astype('datetime64[D]'))

# compute time averages of polygon data
p104_tavg = {}
p001_tavg = {}
for key in p104.keys():
    if not key=='nelement':
        p104_tavg[key] = p104[key].copy()
        p001_tavg[key] = p001[key].copy()
        p104_tavg[key][:] = 0
        p001_tavg[key][:] = 0
for i in range(len(time)):
    #ind = np.logical_and(time>=time_tavg[i], time<(time_tavg[i]+np.timedelta64(1,'D')))
    ind = np.logical_and(time>=(time[i]-np.timedelta64(1,'D')), time<time[i])
    if any(ind):
        for key in p104.keys():
            if not key=='nelement':       
                p104_tavg[key][i] = np.mean(p104[key][ind])
                p001_tavg[key][i] = np.mean(p001[key][ind])

# compute time averages of polygon balance data
p104_bal_tavg = {}
p001_bal_tavg = {}
for key in p104_bal.keys():
    p104_bal_tavg[key] = p104_bal[key].copy()
    p001_bal_tavg[key] = p001_bal[key].copy()
    p104_bal_tavg[key][:] = 0
    p001_bal_tavg[key][:] = 0
    for i in range(len(time_bal)):
        #ind = np.logical_and(time_bal>=(time_bal[i]-np.timedelta64(1,'D')), time_bal<time_bal[i])
        #ind = np.logical_and(time_bal>=(time_bal[i]), time_bal<(time_bal[i]+np.timedelta64(1,'D')))
        ind = np.logical_and(time_bal>(time_bal[i]), time_bal<=(time_bal[i]+np.timedelta64(1,'D')))
        p104_bal_tavg[key][i] = np.mean(p104_bal[key][ind])
        p001_bal_tavg[key][i] = np.mean(p001_bal[key][ind])

# what is the number of elements in the polygons??????????
# in the input file, first polygon is 43 instead of 0
# polygon104 >>> polygon147    has 1701 elements
# polygon1 >>> polygon44       has 951 elements
p001['nelement'] = 950
p104['nelement'] = 1700

# multipliers based on above 
mult_dict = {}
mult_dict['surf'] = 0.1
mult_dict['volume'] = 0.1
mult_dict['totaldepth'] = 0.1
mult_dict['depth'] = 0.1
mult_dict['localdepth'] = 0.1
mult_dict['diat'] = 1
mult_dict['limnutdiat'] = 0.1
mult_dict['limraddiat'] = 0.1
mult_dict['fppdiat'] = 0.1
mult_dict['rcgrodiat'] = 0.1
mult_dict['rcrespdiat'] = 0.1
mult_dict['tavg_dppdiat'] = 1
mult_dict['tavg_rcgrodiat'] = 1
mult_dict['tavg_rcrespdiat'] = 1

# start and end time for plotting
ts = np.datetime64('2013-04-01')
te = np.datetime64('2013-06-01')

sys.exit()


# compare hist-bal and hist measures of productivity
fig, ax = plt.subplots(2,2,figsize=(16,11))
ax1 = ax[0,0]
ax1.set_title('Polygon 104')
ax1.plot(time_bal,p104_bal['Diat,dPPDiat'],label='Diat,dPPDiat')
ax1.plot(time,p104['tavg_dppdiat']*p104_tavg['volume'],label='tavg_dppdiat x volume (using a daily average volume)')
ax1.plot(time_bal+np.timedelta64(1,'D'),p104_bal_tavg['Diat,dPPDiat'],label='Diat,dPPDiat (daily average), time shifted forward one day')
ax1.legend()
ax1.set_xlim((ts,te))
ax1 = ax[1,0]
ax1.plot(time,p104['fppdiat'],'k',label='fppdiat')
ax1.plot(time,p104_tavg['fppdiat'],'b',label='fppdiat (daily rolling average)')
ax1.set_ylabel('fppdiat')
ax2 = ax1.twinx()
ax2.set_ylabel('tavg_fppdiat')
ax2.plot(time,p104['tavg_dppdiat'],'r',label='tavg_dppdiat')
ax1.legend(loc=2)
ax2.legend(loc=3)
ax1.set_xlim((ts,te))
ax1 = ax[0,1]
ax1.set_title('Polygon 1')
ax1.plot(time_bal,p001_bal['Diat,dPPDiat'],label='Diat,dPPDiat')
ax1.plot(time,p001['tavg_dppdiat']*p001_tavg['volume'],label='tavg_dppdiat x volume (using a daily average volume)')
ax1.plot(time_bal+np.timedelta64(1,'D'),p001_bal_tavg['Diat,dPPDiat'],label='Diat,dPPDiat (daily average, time shifted forward one day)')
ax1.legend()
ax1.set_xlim((ts,te))
ax1 = ax[1,1]
ax1.plot(time,p001['fppdiat'],'k',label='fppdiat')
ax1.plot(time,p001_tavg['fppdiat'],'b',label='fppdiat (daily average)')
ax1.set_ylabel('fppdiat')
ax2 = ax1.twinx()
ax2.set_ylabel('tavg_fppdiat')
ax2.plot(time,p001['tavg_dppdiat'],'r',label='tavg_dppdiat')
ax1.legend(loc=2)
ax2.legend(loc=3)
ax1.set_xlim((ts,te))
fig.autofmt_xdate()
fig.suptitle('Compare every measure of productivity for monitoring areas\n(polygons), including both *.his and *-bal.his file data')
plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Compare_Measures_of_Productivity_for_Polygons.png')
plt.close('all')




# compare layer average with layer by layer at stations
for var in var_list:
    fig, ax = plt.subplots(3,2,figsize=(11,11))
    ax1 = ax[0,0]
    ax1.set_title('USGS Station 28')
    ax1.set_ylabel(var)
    ax1.plot(time, np.mean(s028_lay[var], axis=1), label='mean of station 28_i value for i=0-9')
    ax1.legend()
    ax1.set_xlim((ts,te))

    ax1 = ax[1,0]
    ax1.set_ylabel(var)
    ax1.plot(time, s028[var], label='station 28 value')
    ax1.legend()
    ax1.set_xlim((ts,te))

    ax1 = ax[2,0]
    ax1.set_ylabel('navg')
    ratio = np.mean(s028_lay[var], axis=1) / s028[var].values[:,0]
    ax1.plot(time, ratio, '.', label = '(mean of station 28_i) / (station 28 value)')
    if not np.mean(s028[var])==0:
        ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
    ax1.set_xlim((ts,te))

    ax1 = ax[0,1]
    ax1.set_title('USGS Station 35')
    ax1.set_ylabel(var)
    ax1.plot(time, np.mean(s035_lay[var], axis=1), label='mean of station 35_i value for i=0-9')
    ax1.legend()
    ax1.set_xlim((ts,te))

    ax1 = ax[1,1]
    ax1.set_ylabel(var)
    ax1.plot(time, s035[var], label='station 35 value')
    ax1.legend()
    ax1.set_xlim((ts,te))

    ax1 = ax[2,1]
    ax1.set_ylabel('navg')
    ratio = np.mean(s035_lay[var], axis=1) / s035[var].values[:,0]
    ax1.plot(time, ratio, '.', label = '(mean of station 35_i) / (station 35 value)')
    if not np.mean(s028[var])==0:
        ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
    ax1.set_xlim((ts,te))

    fig.autofmt_xdate()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Compare_USGS_Station_Layer_Average_%s.png' % var)
    plt.close('all')



for var in var_list:

    if var in mult_dict.keys():

        if mult_dict[var] == 0.1:
            mult_str = '0.1 x '
        else:
            mult_str = ''
    
        fig, ax = plt.subplots(3,2,figsize=(11,11))
    
        ax1 = ax[0,0]
        ax1.set_title('USGS Station 28 vs. Polygon 104')
        ax1.plot(time, mult_dict[var] * s028[var], label = '%sstation 28 value' % mult_str)
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[1,0]
        ax1.plot(time, p104[var], label = 'polygon 104 value')
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[2,0]
        ratio = p104[var] / (mult_dict[var] * s028[var])
        if np.sqrt(np.var(ratio)) < np.mean(ratio):
            ax1.plot(time, ratio, '.', label = 'polygon 104 value / %sstation 28 value' % mult_str)
            if not np.mean(s028[var])==0:
                ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
            ax1.set_ylabel('ratio (see legend)')
            ax1.legend()
            ax1.set_xlim((ts,te))
        else:
            x = (mult_dict[var] * s028[var])[:,0]
            y = p104[var]
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax1.plot(x, y, '.')  
            ax1.plot(x,p(x), label="y=%.6fx+(%.6f)"%(z[0],z[1]))
            ax1.set_xlabel('%sstation 28 value' % mult_str)
            ax1.set_ylabel('polygon 104 value')  
            ax1.legend()

        ax1 = ax[0,1]
        ax1.set_title('USGS Station 35 vs. Polygon 1')
        ax1.plot(time, mult_dict[var] * s035[var], label = '%sstation 35 value' % mult_str)
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[1,1]
        ax1.plot(time, p001[var], label = 'polygon 1 value')
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
        
        ax1 = ax[2,1]
        ratio = p001[var] / (mult_dict[var] * s035[var])
        if np.sqrt(np.var(ratio)) < np.mean(ratio):
            ax1.plot(time, ratio, '.', label = 'polygon 1 value / %sstation 35 value' % mult_str)
            if not np.mean(s035[var])==0:
                ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
            ax1.set_ylabel('ratio (see legend)')
            ax1.legend()
            ax1.set_xlim((ts,te))
        else:
            x = (mult_dict[var] * s035[var])[:,0]
            y = p001[var]
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax1.plot(x, y, '.')  
            ax1.plot(x,p(x), label="y=%.6fx+(%.6f)"%(z[0],z[1]))
            ax1.set_xlabel('%sstation 35 value' % mult_str)
            ax1.set_ylabel('polygon 1 value')  
            ax1.legend()
    
        fig.autofmt_xdate()
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Compare_USGS_Station_vs_Polygon_%s.png' % var)
        plt.close('all')


# localdepth = depth from water surface to bottom of segment
# totaldepth = total depth water column
# depth = depth of segment

# explore dimensions, depth, surf, volume



# do some special comparisons
fig, ax = plt.subplots(2,2, figsize=(11,8.5))
var = 'limnutdiat'
ax1 = ax[0,0]
ax1.plot(time, p104[var]/p104['nelement'], label = '1/%d x polygon 104 value' % p104['nelement'])
ax1.plot(time, 0.1*s028[var], label = '0.1 x station 28 value')
ax1.legend()
ax1.set_ylabel(var)
ax1.set_title('Polygon 104 has %d elements (per *.inp file)' % p104['nelement'])
ax1 = ax[0,1]
ax1.plot(time, p001[var]/p001['nelement'], label = '1/%d x polygon 1 value' % p001['nelement'])
ax1.plot(time, 0.1*s035[var], label = '0.1 x station 35 value')
ax1.legend()
ax1.set_ylabel(var)
ax1.set_title('Polygon 1 has %d elements (per *.inp file)' % p001['nelement'])
var = 'limraddiat'
ax1 = ax[1,0]
ax1.plot(time, p104[var]/p104['nelement'], label = '1/%d x polygon 104 value' % p104['nelement'])
ax1.plot(time, 0.1*s028[var], label = '0.1 x station 28 value')
ax1.legend()
ax1.set_ylabel(var)
ax1 = ax[1,1]
ax1.plot(time, p001[var]/p001['nelement'], label = '1/%d x polygon 1 value' % p001['nelement'])
ax1.plot(time, 0.1*s035[var], label = '0.1 x station 35 value')
ax1.legend()
ax1.set_ylabel(var)
fig.suptitle('How are limnutdiat and limraddiat aggregated over the monitoring area polygons???')        
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Explore_Polygon_Aggregation_limnutdiat_limraddiat.png')
plt.close('all')

# try normalizing the rates by the number of elements and compare polygons and stations again
for var in ['fppdiat','rcrespdiat','rcgrodiat']:

    if var in mult_dict.keys():

        if mult_dict[var] == 0.1:
            mult_str = '0.1 x '
        else:
            mult_str = ''
    
        fig, ax = plt.subplots(3,2,figsize=(11,11))
    
        ax1 = ax[0,0]
        ax1.set_title('USGS Station 28 vs. Polygon 104')
        ax1.plot(time, mult_dict[var] * s028[var], label = '%sstation 28 value' % mult_str)
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[1,0]
        ax1.plot(time, p104[var]/p104['nelement'], label = '1/%d x polygon 104 value' % p104['nelement'])
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[2,0]
        ratio = p104[var] /p104['nelement'] / (mult_dict[var] * s028[var])
        if np.sqrt(np.var(ratio)) < np.mean(ratio):
            ax1.plot(time, ratio, '.', label = 'polygon 104 value / %sstation 28 value' % mult_str)
            if not np.mean(s028[var])==0:
                ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
            ax1.set_ylabel('ratio (see legend)')
            ax1.legend()
            ax1.set_xlim((ts,te))
        else:
            x = (mult_dict[var] * s028[var])[:,0]
            y = p104[var]/p104['nelement'] 
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax1.plot(x, y, '.')  
            ax1.plot(x,p(x), label="y=%.6fx+(%.6f)"%(z[0],z[1]))
            ax1.set_xlabel('%sstation 28 value' % mult_str)
            ax1.set_ylabel('1/%d x polygon 104 value' % p104['nelement'])  
            ax1.legend()

        ax1 = ax[0,1]
        ax1.set_title('USGS Station 35 vs. Polygon 1')
        ax1.plot(time, mult_dict[var] * s035[var], label = '%sstation 35 value' % mult_str)
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
    
        ax1 = ax[1,1]
        ax1.plot(time, p001[var]/p001['nelement'], label = '1/%d x polygon 1 value' % p001['nelement'])
        ax1.set_ylabel(var)
        ax1.legend()
        ax1.set_xlim((ts,te))
        
        ax1 = ax[2,1]
        ratio = p001[var] /p001['nelement']/ (mult_dict[var] * s035[var])
        if np.sqrt(np.var(ratio)) < np.mean(ratio):
            ax1.plot(time, ratio, '.', label = 'polygon 1 value / %sstation 35 value' % mult_str)
            if not np.mean(s035[var])==0:
                ax1.set_ylim((0.999*np.nanmin(ratio), 1.001*np.nanmax(ratio)))
            ax1.set_ylabel('ratio (see legend)')
            ax1.legend()
            ax1.set_xlim((ts,te))
        else:
            x = (mult_dict[var] * s035[var])[:,0]
            y = p001[var]/p001['nelement'] 
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax1.plot(x, y, '.')  
            ax1.plot(x,p(x), label="y=%.6fx+(%.6f)"%(z[0],z[1]))
            ax1.set_xlabel('%sstation 35 value' % mult_str)
            ax1.set_ylabel('1/%d x polygon 1 value' % p001['nelement'])  
            ax1.legend()
    
        fig.autofmt_xdate()
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Compare_USGS_Station_vs_Polygon_Norm_By_Elements_%s.png' % var)
        plt.close('all')





# revisit hist-bal and try normalizing fppdiat by number of elements
fig, ax = plt.subplots(1,2,figsize=(16,8))
ax1 = ax[0]
ax1.set_title('Polygon 104')
ax1.plot(time_bal,p104_bal['Diat,dPPDiat'],label='Diat,dPPDiat')
ax1.plot(time,p104['tavg_dppdiat']*p104_tavg['volume'],label='tavg_dppdiat x volume (using a daily average volume)')
ax1.plot(time_bal+np.timedelta64(1,'D'),p104_bal_tavg['Diat,dPPDiat'],label='Diat,dPPDiat (daily average, time shifted forward one day)')
ax1.plot(time,p104_tavg['fppdiat']/p104['nelement']*p104['volume'],label='1/%d x fppdiat x volume' % p104['nelement'])
ax1.legend()
ax1.set_xlim((ts,te))
ax1 = ax[1]
ax1.set_title('Polygon 1')
ax1.plot(time_bal,p001_bal['Diat,dPPDiat'],label='Diat,dPPDiat')
ax1.plot(time,p001['tavg_dppdiat']*p001_tavg['volume'],label='tavg_dppdiat x volume (using a daily average volume)')
ax1.plot(time_bal+np.timedelta64(1,'D'),p001_bal_tavg['Diat,dPPDiat'],label='Diat,dPPDiat (daily average, time shifted forward one day)')
ax1.plot(time,p001_tavg['fppdiat']/p001['nelement']*p001['volume'],label='1/%d x fppdiat x volume' % p001['nelement'])
ax1.legend()
ax1.set_xlim((ts,te))
fig.autofmt_xdate()
fig.suptitle('Compare every measure of productivity for monitoring areas\n(polygons), including both *.his and *-bal.his file data')
plt.savefig('../../Plots/troubleshooting/troubleshoot_extra_variable_figs/Compare_Measures_of_Productivity_for_Polygons_NORM_fppdiat.png')
plt.close('all')
