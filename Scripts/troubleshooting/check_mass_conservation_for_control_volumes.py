

########################################################################################
# import python packages
########################################################################################

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import datetime as dt
import geopandas as gpd

########################################################################################
# user input
########################################################################################

# list of balance table csv files, including paths, where you want to check mass balances
balance_table_input = [
'Balance_Tables/diat_Table.csv',
'Balance_Tables/zoopl_v_Table.csv',
'Balance_Tables/zoopl_r_Table.csv',
'Balance_Tables/zoopl_e_Table.csv',
'Balance_Tables/no3_Table.csv',
'Balance_Tables/nh4_Table.csv',
'Balance_Tables/pon1_Table.csv',
'Balance_Tables/detns1_Table.csv',
'Balance_Tables/detps1_Table.csv',
'Balance_Tables/detcs1_Table.csv',
'Balance_Tables/poc1_Table.csv'
]

# folder to save figures
figfolder = 'Balance_Check_CV_Level'

# path to shapefile delineating control volumes
control_volume_shapefile = 'Agg_mod_contiguous_141.shp'

#########################################################################################
# main
#########################################################################################

if not os.path.exists(figfolder):
    os.makedirs(figfolder)

# make lists of keys for flux terms
keys_flux = ['Flux0', 'Flux1', 'Flux2', 'Flux3','Flux4', 'Flux5']

# load just the first data frame in the balance table and count the polygons
df_master = pd.read_csv(balance_table_input[0])
polygons = np.unique(df_master['Control Volume'].values)
npoly = len(polygons)

# loop through the data frames and extract the parameter names into a list
param_list = []
for balance_table in balance_table_input:
    
    # load balance table
    df_master = pd.read_csv(balance_table)
    
    # find the name of the parameter in this balance table (e.g. NO3, Diat) by
    # looking for the first column name that has a comma in it -- the parameter
    # name is to the left of the comma, and append the parameter name to the 
    # list of parameter names
    param = ''
    for cn in df_master.columns:
        if ',' in cn:
            param = cn.split(',')[0]
            break
    assert(not param=='')
    param_list.append(param)
nparams = len(param_list)

# initialize data frames to store errors in transport and dMass/dt for each balance table
df_err_transport = pd.DataFrame(index=np.arange(npoly), columns=param_list, dtype=float)
df_err_transport.index.name = 'polygon'
df_err_masscons = df_err_transport.copy(deep=True)

# loop through the balance tables
for ip in range(nparams):
    
    # get parameter name and name of corresponding balance table
    param = param_list[ip]
    balance_table = balance_table_input[ip]
    
    # load balance table
    df_master = pd.read_csv(balance_table)
    
    # convert time stamp strings to datetime objects, noting that for some reason the format of the 
    # time stamp strings is different in different balance tables -- sometimes Y-M-D and sometimes M/D/Y,
    # so accomodate for this
    try:
        df_master['time'] = [dt.datetime.strptime(timestamp, '%Y-%m-%d').date() for timestamp in df_master['time']]
    except:
        df_master['time'] = [dt.datetime.strptime(timestamp, '%m/%d/%Y').date() for timestamp in df_master['time']]
        
    # extract unique times and find length
    times = np.unique(df_master['time'].values)
    Nt = len(times)
    
    # get reaction keys
    keys_reaction = []
    for col in df_master.columns:
        if (param + ',d') in col:
            keys_reaction.append(col)
    
    # loop through all the polygons
    for polynum in range(npoly):
    
        # select data frame subset for this polygon
        df = df_master[df_master['Control Volume'] == 'polygon%d' % polynum]
        
        # add up the various terms
        Reactions = np.zeros(Nt)
        for key in keys_reaction:
            Reactions += df[key].values.copy()
        Flux_In = np.zeros(Nt)
        for key in keys_flux:
            if any(~np.isnan(df[key].values)):
                Flux_In = Flux_In + df[key].values
        Loads = df[param + ',Loads in'].values
        Transport_In = df[param + ',Transp in'].values
        Transport_Out = df[param + ',Transp out'].values
        Transport = Transport_In + Transport_Out
        dMdt_check = Loads + Transport + Reactions
        if 'Det' in param:
            M = df['Concentration (mg/l)'].values * df['Area'].values * 10.0   # dMass/dt must be calculated differently for detritus
            dMdt = M[1:] - M[0:-1]
            dMdt = np.append(dMdt,dMdt[-1])
        else:
            dMdt = df['dVar/dt'].values
        
        # plot to compare transport in, transport out, and flux out
        plt.figure(figsize=(9,8))
        plt.subplot(3,1,1)
        plt.plot(times,Transport_In/1.0e9,'o',label='Transport In')
        plt.plot(times,Transport_Out/1.0e9,'x',label='Transport Out')
        plt.plot(times,Transport/1.0e9,'s',label='Transport = Transport In + Transport Out')
        plt.plot(times,Flux_In/1.0e9,'--',label='Flux In')
        plt.legend()
        plt.ylabel('Mg/d')
        plt.subplot(3,1,2)
        plt.plot(times,Transport/1.0e9,'o',label='Transport = Transport In + Transport Out')
        plt.plot(times,Flux_In/1.0e9,'--',label='Flux In')
        plt.legend()
        plt.ylabel('Mg/d')
        plt.subplot(3,1,3)
        plt.plot(times,Transport/1.0e9,'o',label='Transport = Transport In + Transport Out')
        plt.plot(times,Reactions/1.0e9,'x',label='Reactions')
        plt.plot(times,Loads/1.0e9,'s',label='Loads')
        plt.plot(times,dMdt/1.0e9,'d',linewidth=2,label='dMass/dt')
        plt.plot(times,dMdt_check/1.0e9,'--',linewidth=2,label='dMass/dt (Check) = Loads + Transport + Reactions')
        plt.legend()
        plt.ylabel('Mg/d')
        plt.suptitle(param + ', Polygon %d' % polynum)
        plt.savefig(figfolder + '/' + param + '_Polygon_%03d.png' % polynum)
        plt.close()
        
        # store errors in transport and in mass conservation in appropriate data frames
        df_err_transport.loc[polynum, param] = np.abs((Flux_In - Transport)).max()/1.0e9
        df_err_masscons.loc[polynum, param] = np.abs((dMdt_check - dMdt)).max()/1.0e9
        
# write errors to data csv files
with open(figfolder + '/_Errors_In_Net_Flux_Compared_to_Net_Transport.csv', 'w') as f:
    f.write('# Max error in Net Flux when compared to Net Transport; units are Mg/d\n')
    df_err_transport.to_csv(f,float_format='%.8f')
with open(figfolder + '/_Errors_In_Mass_Cons_Compared_to_dMdt.csv', 'w') as f:
    f.write('# Max error in Net Loads + Net Transport + Net Reactions when compared to dMdt; units are Mg/d\n')
    df_err_transport.to_csv(f,float_format='%.8f')

# now save images showing maps of transport and mass conservation errors for each parameter
gdf = gpd.read_file(control_volume_shapefile)

# loop through the parameters for the transport errors and mass conservaiton errors, plotting on the grid
for param in param_list:

    # transport errors
    fig, ax = plt.subplots(figsize=(10,10))
    gdf['color'] = df_err_transport.loc[:,param]
    gplot = gdf.plot(ax=ax, column='color', legend=True)
    plt.title(param + ' Error: Net Flux In - Net Transport In (Mg/d)')
    plt.savefig(figfolder + '/_Transport_Error_' + param + '.png')
    plt.close('all')
    
    # mass conservation errors
    fig, ax = plt.subplots(figsize=(10,10))
    gdf['color'] = df_err_masscons.loc[:,param]
    gplot = gdf.plot(ax=ax, column='color', legend=True)
    plt.title(param + ' Error: Net Flux In - Net Transport In (Mg/d)')
    plt.savefig(figfolder + '/_Mass_Conservation_Error_' + param + '.png')
    plt.close('all')

        


