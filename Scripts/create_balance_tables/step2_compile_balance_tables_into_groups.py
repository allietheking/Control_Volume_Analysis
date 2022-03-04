
# building on Zhenlin's control volume analysis
# reads budgets of NH4, NO3, Diatoms, etc. for control volumes and lumps the budgets
# into "groups" A, B, C, etc., computing flux out of each group and reactions within 
# each group. output is a balance table called 'Balance_Table_All_Parameters_By_Group.csv'

#########################################################################################
# user input
#########################################################################################

# directory containting balance tables, including path. note if there are multiple runs,
# there is another layer of run folders inside this directory before the actual balance tables
# water_year = 'WY2018'
run_folders = ['FR18_003']
# output directory, including path
output_dir = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/'

FR = True 


# list of run subfolders inside the directory defined above (if there is only one run and thus no 
# subfolders, set run_folders = [''])
# run_folders = ['cccsd','ebda','ebmud','minus50','plus50','san_jose']
# import os
# run_folders =  [f.name for f in os.scandir(output_dir) if f.is_dir() ]


# names of text files defining groups and their connectivity (includning path)
# group_definition_file   = '/hpcvol1/hpcshared/Scripts/Postprocessors/Control_Volume/control_volume_definitions_141.txt'
# group_connectivity_file = '/hpcvol1/hpcshared/Scripts/Postprocessors/Control_Volume/connectivity_definitions_141.txt'


# for FR Runs (includes new segments ...)
if FR:
    group_definition_file   = '../../Control_Volume_Definitions/control_volume_definitions_FR.txt'
    group_connectivity_file = '../../Control_Volume_Definitions/connectivity_definitions_FR.txt'

# list of balance table csv files, including paths, evaluated at the original control 
# volume level, to read in and compile into groups 
balance_table_list = [
    'diat_Table.csv',
    'zoopl_v_Table.csv',
    'zoopl_r_Table.csv',
    'zoopl_e_Table.csv',
    'no3_Table.csv',
    'nh4_Table.csv',
    'don_Table.csv',
    'pon1_Table.csv',
    'detns1_Table.csv',
    'detns2_Table.csv',
    'oxy_Table.csv',
    'po4_Table.csv',
    'pop1_Table.csv',
    'dop_Table.csv'
    ]

########################################################################################
# import python packages
########################################################################################

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import datetime as dt

###########################################################################################
# load dictionaries that define the groups (A, B, C, etc.) by the list of control volumes
# that make up each group (group_dict) and the list of pairs of control volumes making up
# fluxes out of the N, S, E, and W faces of each group (flux_pairs_dict)
###########################################################################################

# read from text file a dictionary grouping control volumes into lettered groups 
# format:
#   GROUP_NAME : CV1, CV2, CV3, ...
# where:
#   GROUP_NAME = name of group
#   CV1, CV2, CV3, ... = numbers of control volumes comprising the group
group_dict = {}
with open(group_definition_file,'r') as f: 
    for line in f.readlines():
        # split line at ':' into key and control volume list, dropping newline character first
        (key, val) = line.rstrip('\n').split(':')  
        # remove whitespace surrounding key 
        key = key.strip()
        # add list to dictionary
        print(key)
        group_dict[key] = eval(val) 
      
        
# read from text file a dictionary specifying all the fluxes between 
# control volumes making up fluxes from each lettered group to the north, 
# south, east, and west
# format:
#   LETTER2DIRECTION : [[CVf,CVt]_1, [CVf,CVt]_2, [CVf,CVt]_3, ...]
# where:
#   LETTER = name of group
#   DIRECTION = N, S, E, or W stand for "to the north", south, east, or west, respectively
#   [CVf, CVt]_i = defines a flux from control volume CVf to control volume CVt
#   where i = 1, 2, 3, ... make up all the fluxes for a particular group/direction pair
flux_pairs_dict = {}
with open(group_connectivity_file,'r') as f: 
    for line in f.readlines():
        # split line at ':' into key and connectivity volume list, dropping newline character first
        (key, val) = line.rstrip('\n').split(':') 
        # trim whitespace around key
        key = key.strip()
        # add list to dictionary
        flux_pairs_dict[key] = eval(val) 
        
# loop through the run folders 
for run_folder in run_folders:
    print(run_folder)
    #######################################################################################
    # Read control volume budgets from the balance tables, extract time info, 
    # and do some basic checks to make sure the time stamps and control volumes 
    # line up so we can combine this data into one big data frame  
    ########################################################################################
    
    # make list of paths to balance tables in this run folder
    balance_table_input = [os.path.join(output_dir,run_folder,bt) for bt in balance_table_list]
    
    # read in each balance table to a data frame, and compile all these data frames into a  
    # list of data frames called df_params
    Nparams = len(balance_table_input)
    df_params = []
    for balance_table_file in balance_table_input:
        df = pd.read_csv(balance_table_file)
        df_params.append(df)
    
    # convert time stamp strings to datetime objects, noting that for some reason the format of the 
    # time stamp strings is different in different balance tables -- sometimes Y-M-D and sometimes M/D/Y,
    # so accomodate for this
    for df in df_params:
        try:
            df['time'] = [dt.datetime.strptime(timestamp, '%Y-%m-%d').date() for timestamp in df['time']]
        except:
            df['time'] = [dt.datetime.strptime(timestamp, '%m/%d/%Y').date() for timestamp in df['time']]
    
    # extract unique dates and note number of time steps, and also create a date string
    # corresponding to these dates
    dates = np.sort(np.unique(df_params[0]['time'].values))
    Nt = len(dates)
    datestrings = [dt.datetime.strftime(date,'%Y-%m-%d') for date in dates]
    
    # the tables should all be the same size, and the time stamps and polygons should match -- check
    # that this is the case, and if not, throw an error
    for ip in np.arange(1,Nparams):
        assert(all(df_params[0]['Control Volume']==df_params[ip]['Control Volume']))
        assert(all(df_params[0]['time']==df_params[ip]['time']))
    
    ##########################################################################################
    # Intialize a gigantic master data frame to store all of the data by group (A, B, C, ...)
    ##########################################################################################
    
    # first note how many groups there are, and make sure we have the right number of connections 
    # defined in our flux pair dictionary -- N/S/E/W (4 connections) for each group
    Ng = len(group_dict)
    assert(Ng*4 == len(flux_pairs_dict))
    
    # initialize a list that will contain all the column names for the master data frame, 
    # starting with the variables that are the same across all the parameters
    column_names = ['group','time','Volume (m^3)','Area (m^2)']
    
    # compile a list of parameter names, one corresponding to each balance table
    param_name_list = []
    
    # compile a list of the reaction term names for each parameter -- this will be a list 
    # of lists
    param_reaction_name_list = []
    
    # to get the rest of the column names, need to look at the names of the columns in each
    # of the balance table data frames
    for df in df_params:
        
        # find the name of the parameter in this balance table (e.g. NO3, Diat) by
        # looking for the first column name that has a comma in it -- the parameter
        # name is to the left of the comma, and append the parameter name to the 
        # list of parameter names
        param_name = ''
        for cn in df.columns:
            if ',' in cn:
                param_name = cn.split(',')[0]
                param_name_list.append(param_name)
                break
        assert(not param_name=='')
        
        # add column names that are the same for each parameter and don't necessarily
        # correspond to column names in the original balance table csv files
        column_names.extend([param_name + ',' + 'Mass (Mg)', 
                             param_name + ',' + 'dMass/dt (Mg/d)',
                             param_name + ',' + 'Net Flux In (Mg/d)',
                             param_name + ',' + 'Net Load (Mg/d)',
                             param_name + ',' + 'Net Transport In (Mg/d)',
                             param_name + ',' + 'Net Reaction (Mg/d)',
                             param_name + ',' + 'dMass/dt, Balance Check (Mg/d)',
                             param_name + ',' + 'Flux In from N (Mg/d)',
                             param_name + ',' + 'Flux In from S (Mg/d)',
                             param_name + ',' + 'Flux In from E (Mg/d)',
                             param_name + ',' + 'Flux In from W (Mg/d)'])
                
        # note number of columns in the data frame for this parameter
        Nc = len(df.columns)
                             
        # append the reaction terms, which are unique to each parameter, also save the 
        # list of reaction names for each parameter. this code assumes that the reaction
        # terms are all the column names with the parameter name followed by a comma and 
        # the letter "d":
        reaction_name_list = []
        for rt in df.columns:
            if (param_name + ',d') in rt:
                column_names.append(rt + ' (Mg/d)')
                reaction_name_list.append(rt + ' (Mg/d)')
        param_reaction_name_list.append(reaction_name_list)
    
    # now that we have all the column names compiled, initialize the master data frame
    df_master = pd.DataFrame(columns=column_names)
    
    # fill in the groups and timestamps in the master data frame such that the order
    # of the groups are A, A, A, ...., B, B, B, ...., C, C, C, ...., etc. with the number 
    # of repeats of each gropup letter equal to the number of time stamps
    group_list = []
    datestrings_list = []
    for group in group_dict.keys():
        group_list += [group] * Nt            
        datestrings_list += datestrings.copy() 
    df_master['group'] = group_list.copy()
    df_master['time'] = datestrings_list.copy()
    
    # initalize all of the other columns in the data frame, i.e., everything besides
    # group and time stamp, as arrays of zeros of type float
    for col in range(2,len(df_master.columns)):
        df_master.iloc[:,col]=np.zeros(Nt*Ng, dtype=float)
       
    ########################################################################################
    # for each group/face pair, i.e., A2N, A2S, ... , L2E, L2W, add up all of the
    # fluxes between pairs of control volumes that make up the total flux through 
    # that face
    ########################################################################################
    
    # loop through all of the group/face pairs (A2N, A2S, ... , L2E, L2W)
    for flux_key in flux_pairs_dict.keys():    
    
        # extract group letter (A, B, C, ... ) and face (N, S, E, W) from flux_key
        group = str(flux_key).split()[0]
        NSEW = str(flux_key).split()[-1]
        
        # initialize total flux through this group/face pair at zero for each of the parameters,
        # storing the fluxes in a list
        flux_params = []
        for ip in range(Nparams):
            flux_params.append(np.zeros(Nt, dtype=float))
        
        # loop through pairs of control volumes comprising total flux through this group/face, 
        # and add up the fluxes
        for flux_pair in flux_pairs_dict[flux_key]:  
            
            # each flux_pair is a 2-item list where the first item is the "from" 
            # control volume/polygon and the second item is the "to" control volume/polygon
            print(flux_pair)
            from_poly = flux_pair[0]
            to_poly = flux_pair[1]        
            
            # using logical indexing, extract the rows of the balance table data frames corresponding to 
            # the "from" polygon in this particular flux pair for each parameter, storing the values in a list
            rows_from_param = []
            for ip in range(Nparams):
                rows_from_param.append(df_params[ip].loc[df_params[ip]['Control Volume'] == ('polygon%d' % from_poly)])
            
            # only need one row to identify which control volume is the "to" polygon, so pick out the first row
            row_to_param = []
            for ip in range(Nparams):
                row_to_param.append(rows_from_param[ip].iloc[0])
                
            # count the number of fields labeled "To_polyN" where N is an integer
            n2poly = 0
            for key in row_to_param[0].keys():
                if 'To_poly' in key:
                    n2poly = n2poly + 1
            
            # note that the way Zhenlin handled the "to" control volumes is by creating eight fields called 
            # "To_polyX" where X = 0, 1, ..., 7, that map to the numbers of all the "to" control volumes for each
            # "from" polygon. using our one example row, find the number X of the "to" polygon in the corresponding 
            # to the "to" control volume in our flux pair -- should be the same for all substances but doing them 
            # separately to narrow down origin of any potential errors
            polynum_param = []
            for ip in range(Nparams):
                polynum_param.append(None)
            for p in range(0,n2poly):
                for ip in range(Nparams):
                    if row_to_param[ip]['To_poly%d' % p] == to_poly:
                        polynum_param[ip] = p
            
            # if we can't find the "to" polygon, polynum will remain None, and the 
            # following will throw an assertion error. if you get this error, it means
            # there is either an error in the definition of the flux dictionary or an error
            # in Zhenlin's data tables, probably the dictionary
            for ip in range(Nparams):
                assert(polynum_param[ip]+1)
            
            # add the flux for this control volume pair flux to the total flux for the group/side pair
            for ip in range(Nparams):
                flux_params[ip] += rows_from_param[ip]["Flux%d" % polynum_param[ip]].values.copy()
            
        # add the total flux for this group/side pair to the master data frame in the appropriate place
        # switching direction because while Zhenlin's "To_poly" implies these fluxes are out, they are
        # actually in. divide by 1e6 to converg g/d to Mg/d
        ind_m = df_master['group'] == group
        for ip in range(Nparams):
            df_master.loc[ind_m,param_name_list[ip] + ',Flux In from %s (Mg/d)' % NSEW] = flux_params[ip].copy()/1.0e6
           
    ########################################################################################
    # For each group (A, B, C, ...) add up all the contributions from the multiple
    # control volumes that comprise the group to get the total area, volume, rate of 
    # change of mass, loads, and reaction rates
    ########################################################################################
    
    # loop through the groups
    for group in group_dict.keys():
    
        # loop through all the parameters
        for ip in range(Nparams):
            
            # get the dataframe corresponding to this parameter and get the parameter name
            df = df_params[ip]
            param_name = param_name_list[ip]
            reaction_name_list = param_reaction_name_list[ip]
    
            # if this is the first parameter, initialize the area and volume of the group, 
            # which are the same for all the parameters, at zero
            if ip==0:
                Volume = np.zeros(Nt, dtype=float)
                Area = np.zeros(Nt, dtype=float)
        
            # initialize the variables that are the same for all parameters at zero
            Mass = np.zeros(Nt, dtype=float)
            dMdt = np.zeros(Nt, dtype=float)
            Load_In = np.zeros(Nt, dtype=float)
            Transport_In = np.zeros(Nt, dtype=float)
            
            # initialize all the reaction terms at zero
            Nr = len(reaction_name_list)
            reaction_term_values = []
            for ir in range(Nr):    
                reaction_term_values.append(np.zeros(Nt, dtype=float))
            
            # loop through all the control volumes that comprise this group, and add 
            # up the reaction terms, masses, volumes, loads, etc.
            for cvnum in group_dict[group]:
            
                # logical index for selecting rows of substance-specific data frames 
                # corresponding to this control volume (a.k.a. polygon) 
                ind_s = df['Control Volume'] == 'polygon%d' % cvnum
               
                # if this is the first parameter, add volumes and areas
                if ip==0:
                    Volume    += df.loc[ind_s,'Volume'].values.copy()
                    Area      += df.loc[ind_s,'Area'].values.copy()
                
                # add the mass, multiplying concentration by volume to get mass
                Mass  += df.loc[ind_s,'Volume'].values.copy() * df.loc[ind_s,'Concentration (mg/l)'].values.copy()
            
                # add the mass rate of change, loads, and transport in
                dMdt         +=   df.loc[ind_s,'dVar/dt'].values.copy()
                Load_In      += ( df.loc[ind_s,param_name + ',Loads in'].values.copy() + df.loc[ind_s,param_name + ',Loads out'].values.copy() )
                Transport_In += ( df.loc[ind_s,param_name + ',Transp in'].values.copy() + df.loc[ind_s,param_name + ',Transp out'].values.copy() )
                
                # add the reaction rates, noting that in dataframe df, the reaction name doesn't
                # include the units so need to trim the reaction name string
                for ir in range(Nr):
                    reaction_term_values[ir] += df.loc[ind_s,reaction_name_list[ir][0:-7]].values.copy()
                
            # assign totals to appropriate group in master data frame
            ind_m = df_master['group'] == group
            
            # if this is the first paramter, add volumes and areas
            if ip==0:
                df_master.loc[ind_m,'Volume (m^3)']  = Volume.copy()
                df_master.loc[ind_m,'Area (m^2)']    = Area.copy()
            
            # set mass, mass rate of change, net load, net transport in, dividing by 1.0e6 to 
            # convert g to Mg and g/d to Mg/d
            df_master.loc[ind_m,param_name + ',Mass (Mg)']               = Mass.copy()/1.0e6
            df_master.loc[ind_m,param_name + ',dMass/dt (Mg/d)']         = dMdt.copy()/1.0e6
            df_master.loc[ind_m,param_name + ',Net Load (Mg/d)']         = Load_In.copy()/1.0e6
            df_master.loc[ind_m,param_name + ',Net Transport In (Mg/d)'] = Transport_In.copy()/1.0e6
            
            # set reaction terms, dividing by 1.0e6 to convert g to Mg and g/d to Mg/d
            for ir in range(Nr):
                df_master.loc[ind_m,reaction_name_list[ir]] = reaction_term_values[ir].copy()/1.0e6
    
    
    ########################################################################################
    # Add up the fluxes and reactions for each parameter to get the net and add to dataframe
    ########################################################################################
    
    # loop through the parameters
    for ip in range(Nparams):
    
        # get the name of the parameter and the list of reaction names for this parameter
        param_name = param_name_list[ip]
        reaction_name_list = param_reaction_name_list[ip]
        
        # add N,S,E,W components to get the net flux out and set value in the master data
        # frame
        df_master[param_name + ',Net Flux In (Mg/d)'] = (  
                                  df_master[param_name + ',Flux In from N (Mg/d)'].values.copy()
                                + df_master[param_name + ',Flux In from S (Mg/d)'].values.copy()
                                + df_master[param_name + ',Flux In from E (Mg/d)'].values.copy()
                                + df_master[param_name + ',Flux In from W (Mg/d)'].values.copy() )
        
        # add up all the reaction rates to get the net reaction rate and set value in 
        # the master data frame
        Nr = len(reaction_name_list)
        net_reaction = df_master[reaction_name_list[0]].values.copy()
        for ir in range(1,Nr):
            net_reaction += df_master[reaction_name_list[ir]].values.copy()
        df_master[param_name + ',Net Reaction (Mg/d)'] = net_reaction.copy()
        
        # do a balance check on the rate of change of mass
        df_master[param_name + ',dMass/dt, Balance Check (Mg/d)'] = ( 
                                  df_master[param_name + ',Net Load (Mg/d)'].values.copy()
                                + df_master[param_name + ',Net Reaction (Mg/d)'].values.copy()  
                                + df_master[param_name + ',Net Transport In (Mg/d)'].values.copy() )
    
        
    ########################################################################################
    # Print master dataframe to *.csv
    ########################################################################################
    
    df_master.to_csv(os.path.join(output_dir,run_folder,
                                   'Balance_Table_All_Parameters_By_Group.csv'),index=False)
        
        
    



