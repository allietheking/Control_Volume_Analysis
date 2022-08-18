
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

# to test mass conservation, let's just look at polygon 22, which is bordered by 
# transects 84, 58, 81, 83 (and possibly also small transects 80, 82, 97, 95)
test_polygon = 22
test_transects = [58, 80, 81, 82, 83, 84]
test_transects_direction = [1, -1, -1, -1, -1, -1]   # 1 for clockwise, -1 for counter clockwise

# run folder 
#run_folder = 'G141_13to18_182'
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

# read the shapefiles and isolate the test polygon and adjacent transects
gdf     = gpd.read_file(shpfn_tran).iloc[test_transects]
poly_df = gpd.read_file(shpfn_poly).iloc[test_polygon:test_polygon+1]

# make a plot to make sense of the test polygon and adjacent transect directions
for i in range(len(gdf)):

    geo = gdf.iloc[i].geometry
    xy = geo.coords.xy
    x1 = xy[0][0]
    y1 = xy[1][0]
    
    if test_transects_direction[i]==1:
        dirstr = 'clockwise'
    else:
        dirstr = 'COUNTERclockwise'

    fig, ax = plt.subplots()
    poly_df.plot(ax=ax,edgecolor='k')
    gdf.iloc[i:i+1].plot(ax=ax,color='r')
    ax.plot(x1,y1,'ro')
    ax.set_title('transect%d, first point is red dot\ndirection = %s' % (test_transects[i],dirstr))

plt.close('all')

# if delta, add "delta" prefix to the output folder
run_folder_out = run_folder
if is_delta:
    run_folder_out = 'Delta_' + run_folder_out
    
# path to his and his bal files
histfn      = os.path.join(path,run_folder,'dwaq_hist.nc')
histbal_fn  = os.path.join(path,run_folder,'dwaq_hist_bal.nc')

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

# loop through all the parameters (nh4, no3, diat, etc.)
varnames = [var.lower() for var in hbdata.sub.values]
for varname in ['NH4']:#varnames:

    
    # transect name list
    transect_name_list = ['transect%04d' % t for t in test_transects]

    ##%% Get the data for all the transects 
    TransectBL = [(name in transect_name_list) for name in hdata.location_names.values[0]]
    indT = np.where(TransectBL)[0]
    varT = hdata.isel(nSegment=indT)[varname.lower()]
    
    # get the his data just for the test polygon
    PolygonBL = [(name=='polygon%d' % test_polygon) for name in hdata.location_names.values[0]]
    indP = np.where(PolygonBL)[0]
    varP = hdata.isel(nSegment=indP)[varname.lower()]
    Vp = hdata.isel(nSegment=indP)['volume']  
    
    # get the bal data just for the test polygon
    PolygonBL_bal = [(name=='polygon%d' % test_polygon) for name in hbdata.region.values]
    indP_bal = np.where(PolygonBL_bal)[0]
    fieldBL = [varname.lower()+',' in name.lower() for name in hbdata.field.values]
    indF = np.where(fieldBL)[0]
    varP_bal = hbdata.isel(region=indP_bal).isel(field=indF)

    # get the time axes for the balance and hist files
    time_his = hdata.time.values
    time_bal = hbdata.time.values

    # compute a conversion factor to change units to per day
    cf_his = np.timedelta64(1,'D')/(time_his[1]-time_his[0])
    cf_bal = np.timedelta64(1,'D')/(time_bal[1]-time_his[0])

    # give the polygon concentraiton and volume friendlier names and friendlier dimensions, then compute mass
    conc = varP.values[:,0]
    volume = Vp.values[:,0]
    mass = conc*volume

    # compute the flux in across the transects, converting to 1/d
    flux = varT.values[:,0]*test_transects_direction[0] * cf_his
    for i in range(1,len(test_transects)):
        flux = flux + varT.values[:,i]*test_transects_direction[i] * cf_his

    # compute the net load and the net transport, and net reaction from the balance, converting to 1/d
    loads = (varP_bal.sel(field='NH4,Loads in').bal.values[:,0] + varP_bal.sel(field='NH4,Loads out').bal.values[:,0]) * cf_bal
    transport = (varP_bal.sel(field='NH4,Transp in').bal.values[:,0] + varP_bal.sel(field='NH4,Transp out').bal.values[:,0]) * cf_bal

    # get a list of reactions from the balance file
    rx_list = []
    for field in varP_bal.field.values:
        if ',d' in field:
            rx_list.append(field)

    # compute net reaction from balance file, converting to 1/d
    reaction = varP_bal.sel(field=rx_list[0]).bal.values[:,0] * cf_bal
    for rx in rx_list[1:]:
        reaction = reaction + varP_bal.sel(field=rx).bal.values[:,0] * cf_bal

    # experiment with time averaging windows to get the net flux in across the transects to line up with the net transport in
    # from the balance file ... this one is the winner:
    flux_bal = np.zeros(np.shape(time_bal))
    for i in range(len(time_bal)):
        ind = np.logical_and(time_his>(time_bal[i]-np.timedelta64(6,'h')), time_his<=time_bal[i])
        flux_bal[i] = np.mean(flux[ind])

    # plot to compare net flux across transects with net transport from balance output
    fig, ax = plt.subplots()
    ax.plot(time_his, flux, label='net flux across transects from *.his file')
    ax.plot(time_bal, transport, label='net transport in from *-his.bal file')
    ax.plot(time_bal, flux_bal, '--', label='time average net flux, using t>t-dt, t<=t (backwards average, closed on the right)')
    ax.legend()

    # compute dM/dt from the mass balance
    dmdt_bal = loads + transport + reaction

    # try to compute equivalent change in mass from *.his output
    dmdt_his = np.zeros(len(time_bal))
    for i in range(1,len(time_bal)):
        ind1 = time_his==time_bal[i-1]
        ind2 = time_his==time_bal[i]
        dmdt_his[i] = (mass[ind2] - mass[ind1]) * cf_bal

    # plot to compare dm/dt
    fig, ax = plt.subplots()
    ax.plot(time_bal, dmdt_bal, label='dM/dt that closes terms in *-his.bal file')
    ax.plot(time_bal, dmdt_his, label='dM/dt from *.his file, using M(t) - M(t-dt)')

    # now try resampling onto a daily axis, see if we can do it correctly...
    varP1 = varP.resample(time='1D').nearest()
    Vp1 = Vp.resample(time='1D').nearest()
    Vp_mean1 = Vp.resample(time='1D',closed='right',label='right').mean(dim='time')
    varT1 = varT.resample(time='1D',closed='right',label='right').sum(dim='time')
    varP_bal1 = varP_bal.resample(time='1D',closed='right',label='right').sum(dim='time')

    # check if the resampled data conserves mass
    time_his1 = varP1.time.values
    time_bal1 = varP_bal1.time.values 
    
    # give the polygon concentraiton and volume friendlier names and friendlier dimensions, then compute mass
    conc1 = varP1.values[:,0]
    volume1 = Vp1.values[:,0]
    mass1 = conc1*volume1

    # compute the flux in across the transects
    flux1 = varT1.values[:,0]*test_transects_direction[0] 
    for i in range(1,len(test_transects)):
        flux1 = flux1 + varT1.values[:,i]*test_transects_direction[i]

    # compute the net load and the net transport, and net reaction from the balance
    loads1 = (varP_bal1.sel(field='NH4,Loads in').bal.values[:,0] + varP_bal1.sel(field='NH4,Loads out').bal.values[:,0])
    transport1 = (varP_bal1.sel(field='NH4,Transp in').bal.values[:,0] + varP_bal1.sel(field='NH4,Transp out').bal.values[:,0]) 

    # compute net reaction from balance file, converting to 1/d
    reaction1 = varP_bal1.sel(field=rx_list[0]).bal.values[:,0] 
    for rx in rx_list[1:]:
        reaction1 = reaction1 + varP_bal1.sel(field=rx).bal.values[:,0] 

    # compute dM/dt from the mass balance
    dmdt_bal1 = loads1 + transport1 + reaction1

    # compute dM/dt from the hist data (note it is already per day), and add a zero in front
    dmdt_his1 = mass1[1:] - mass1[0:-1]
    dmdt_his1 = np.insert(dmdt_his1, 0, 0)


 

    sys.exit()

    

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
