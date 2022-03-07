
'''
This script converts *.his and *-bal.his data to *.csv formatted balance tables containing
daily fluxes and reaction terms in the "monitoring regions" defined by polygons and transects.
Note this script requires a python environment with geopandas installed. On HPC load Zhenlin's 
environment by executing the following from the command line: 
    source activate /home/zhenlin/.conda/envs/my_root
'''

################
### USER INPUT
################

run_folders = ['G141_13to18_88']
water_year = 'WY2018'
FR = False


# output directory, including path
output_dir = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/'


# path to folder that contains run folders, which in turn contain model output including dwaq_hist.nc 
# and *-bal.his files. if there are no run subfolders, this is the direct path to dwaq_hist.nc and *-bal.his + 
# path to shapefiles defining polygons and transects for the monitoring regions (make sure these are the same ones used to
# initialize the DWAQ run)
if FR:
    shpfn_tran =  '/richmondvol1/hpcshared/inputs/shapefiles/Agg_exchange_lines_plus_subembayments_shoal_channel.shp' # used in FR17003, FR13003
    shpfn_poly = '/richmondvol1/hpcshared/inputs/shapefiles/Agg_mod_contiguous_plus_subembayments_shoal_channel.shp'
    path = '/richmondvol1/hpcshared/Full_res/%s//' % water_year
    print('FULL RESOLUTION SHAPEFILES.')
else:
    shpfn_poly =  '/richmondvol1/hpcshared/inputs/shapefiles/Agg_mod_contiguous_141.shp'
    shpfn_tran = '/richmondvol1/hpcshared/inputs/shapefiles/Agg_exchange_lines_141.shp'
    path = '/richmondvol1/hpcshared/Grid141/WY13to18//'
    print('AGGREGATED SHAPEFILES.')







##################
# IMPORT MODULES
##################

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from matplotlib.ticker import FormatStrFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import pandas as pd
from matplotlib.patches import Rectangle
sys.path.append('/richmondvol1/siennaw/lib/stompy//')
import stompy.model.delft.io as dio
import datetime 
try:
    import geopandas as gpd
except:
    raise Exception('\n' + 
          'ERROR: this Python environment does not include geopandas.\n' + 
          'To load environment with geopandas, exit Python, execute the\n' + 
          'following at the command line:\n' + 
          '    source activate /home/zhenlin/.conda/envs/my_root\n' + 
          'and then re-launch Python.\n')
    
##################
# FUNCTIONS
##################


user = os.getlogin()
scriptname= __file__ 
today= datetime.datetime.now().strftime('%b %d, %Y')
READ_ME_txt = 'README ! \nThis script was produced on %s by %s using this script:\n %s\nFolder: %s' % (today, user, scriptname, os.getcwd())

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


# richmondvol1/hpcshared/Grid141/WY13to18//G141_13to18_74/dwaq_hist.nc
# h:\hpcshared\Grid141\WY13to18\G141_13to18_74\dwaq_hist.nc

# loop through run folders
for run_folder in run_folders:


    
    # path to his and his bal files
    histfn      = os.path.join(path,run_folder,'dwaq_hist.nc')
    histbal_fn  = os.path.join(path,run_folder,'dwaq_hist_bal.nc')
    print(histfn)
    # open hist file
    # try:
    hdata = xr.open_dataset(histfn)
    
    # except:
    #     print('error')
    #     print(histfn)
    #     continue 
    
    # open hist-bal file
    try:
        hbdata = xr.open_dataset(histbal_fn)
    except:
        fn = os.path.join(path,run_folder,'sfbay_dynamo000-bal.his')
        hbdata = dio.bal_his_file_xarray(fn)
        print('it worked ... ') 
        hbdata.to_netcdf(histbal_fn)
        print('it worked ... 1') 
        hbdata = xr.open_dataset(histbal_fn)
        
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
    print(hbdata.region.values)
    if polyc < 2:
        print('skipping run : %s' % histbal_fn)
        continue 
    
    print('Reading in %s' % run_folder)
    # create the balance tables directory if it doesn't exist yet
    if not os.path.exists(os.path.join(output_dir,run_folder)):
        os.makedirs(os.path.join(output_dir,run_folder))

        
    # loop through all the parameters (nh4, no3, diat, etc.)
    varnames = [var.lower() for var in hbdata.sub.values]
    for varname in varnames:
    
        # create name for balance table output file
        outfile = os.path.join(output_dir,run_folder,varname +'_Table.csv')
        outREADME = os.path.join(output_dir,run_folder,'README.txt')
        
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
        
        # check the frequency of the data, and if it is hourly, resample it onto a daily axis
        deltat_P = (varP.time[1]-varP.time[0]).values
        deltat_T = (varT.time[1]-varT.time[0]).values
        deltat_B = (varP_bal.time[1]-varP_bal.time[0]).values
        # ... polygon output is instantaneous
        if deltat_P==np.timedelta64(3600000000000):
            varP = varP.resample(time='D').nearest()
        elif deltat_P==np.timedelta64(86400000000000):
            pass
        else:
            raise Exception('ERROR: his file must have daily or hourly output time stamp')
        # ... transect output should be integrated in time, bizarrely the integral is open on the 
        # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
        # ... so add a day to it
        if deltat_T==np.timedelta64(3600000000000):
            varT = varT.resample(time='D',closed='right').sum(axis=0)[0:-1,:]
            varT['time'] = varT['time'] + np.timedelta64(86400000000000,'ns')
        elif deltat_T==np.timedelta64(86400000000000):
            pass
        else:
            raise Exception('ERROR: his file must have daily or hourly output time stamp')
        # ... balance output should be integrated in time, bizarrely the integral is open on the 
        # ... left and closed on the right, i.e., integral is from 0<t<=T. note the time ends up shifted
        # ... so add a day to it
        if deltat_B==np.timedelta64(3600000000000):
            varP_bal = varP_bal.resample(time='D').sum(axis=0)[0:-1,:]
            varP_bal['time'] = varP_bal['time'] + np.timedelta64(86400000000000,'ns')
        elif deltat_B==np.timedelta64(86400000000000):
            pass
        else:
            raise Exception('ERROR: his-bal file must have daily or hourly output time stamp')
        
        # now make sure polygon, transect, and balance data have same number of time steps (this 
        # condition may be violated in case of incomplete simulation)
        print(varP.time.values[0])
        print('^^')
        tmin = np.min([varP.time.values[-1],varT.time.values[-1],varP_bal.time.values[-1]])
        varP = varP.where(varP.time<=tmin,drop=True)
        varT = varT.where(varT.time<=tmin,drop=True)
        varP_bal = varP_bal.where(varP_bal.time<=tmin,drop=True)
        
        # load shapefiles
        gdf     = gpd.read_file(shpfn_tran)
        left    = gdf.left
        right   = gdf.right
        pmax = max(left.max(),right.max()) +1 # total number of polygons and they are zero-based
        
        poly_df = gpd.read_file(shpfn_poly)
        
        Vp = hdata.isel(nSegment=indP)['volume']  
        Vp = Vp.resample(freq='D', dim = 'time', how='mean', skipna = True)
        # Vp_mean = np.array(Vp.mean(axis=0))
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
                column_list = ['Control Volume','Concentration (mg/l)','Volume',
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
            pi_df_comb['Volume'] = Vp.values[:,i] #Vp_mean[i]
            pi_df_comb['Area'] = Area[i]    
            
            
            df_output = pd.concat([df_output,pi_df_comb])
                
        df_output = df_output[column_list]    
        
        #%% 
            
        df_output.to_csv(outfile,columns=column_list)   
        print('Saved %s' % outfile)
        
        with open(outREADME, 'w') as f:
            f.write(READ_ME_txt)
         

