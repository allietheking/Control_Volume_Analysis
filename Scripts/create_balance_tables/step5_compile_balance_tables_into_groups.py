
# building on Zhenlin's control volume analysis
# reads budgets of NH4, NO3, Diatoms, etc. for control volumes and lumps the budgets
# into "groups" A, B, C, etc., computing flux out of each group and reactions within 
# each group. output is a balance table called 'Balance_Table_All_Parameters_By_Group.csv'

#########################################################################################
# user input
#########################################################################################

# directory containting balance tables, including path. note if there are multiple runs,
# there is another layer of run folders inside this directory before the actual balance tables
#run_folders = ['G141_13to18_125'] 
run_folders = ['FR13_003','FR13_023','FR17_003','FR17_017','FR18_005'] 

# output directory, including path
output_dir = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/'

# is full resolution?
FR=True

# list of run subfolders inside the directory defined above (if there is only one run and thus no 
# subfolders, set run_folders = [''])
# run_folders = ['cccsd','ebda','ebmud','minus50','plus50','san_jose']
# import os
# run_folders =  [f.name for f in os.scandir(output_dir) if f.is_dir() ]


# names of text files defining groups and their connectivity (includning path)
# group_definition_file   = '/hpcvol1/hpcshared/Scripts/Postprocessors/Control_Volume/control_volume_definitions_141.txt'
# group_connectivity_file = '/hpcvol1/hpcshared/Scripts/Postprocessors/Control_Volume/connectivity_definitions_141.txt'


# for FR Runs (includes new segments defined by sienna and whole bay added by allie...)
if FR:
    group_definition_file   = '../../Definitions/control_volume_definitions_FR.txt'
    group_connectivity_file = '../../Definitions/connectivity_definitions_FR.txt'
# for AGG runs (includs now whole bay group added by allie)
else:
    group_definition_file   = '../../Definitions/control_volume_definitions_141.txt'
    group_connectivity_file = '../../Definitions/connectivity_definitions_141.txt'

# list of parameters corresponding to balance tables
param_list = ['DIN','Algae','Diat','Green','NO3','NH4','TN','TN_include_sediment']
Nparams = len(param_list)

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

# note how many groups there are, and make sure groups are defined only once, also make
# sure there are no repeat connections 
Ng = len(group_dict)
assert Ng == len(np.unique(list(group_dict.keys())))
assert len(list(flux_pairs_dict.keys())) == len(np.unique(list(flux_pairs_dict.keys())))

# find the number of unique types of connections (N/S/E/W also possibly NE/SW/etc)
direction_list = []
for key in flux_pairs_dict.keys():
    direction = key[key.find(' to ')+4:]
    if not direction in direction_list:
        direction_list.append(direction)
print('flux directions include: ')
for direction in direction_list:
    print('    ' + direction)
        
# loop through the run folders 
for run_folder in run_folders:

    print(run_folder)

    # loop through the parameters
    for ip in range(Nparams):

        # get the current parameter
        param_name = param_list[ip]

        #######################################################################################
        # Read control volume budgets from the balance table and extract time info
        ########################################################################################

        # read in the balance table corresponding to this parameter
        balance_table_fn = param_name.lower() + '_Table.csv'
        try:
            df = pd.read_csv(os.path.join(output_dir,run_folder,balance_table_fn))
        except:
            print('error reading %s, skipping this one ...' % balance_table_fn)
            continue
        
        # print
        print('Building grouped balance table for %s ...' % param_name)

        # convert time stamp strings to datetime objects, noting that for some reason the format of the 
        # time stamp strings is different in different balance tables -- sometimes Y-M-D and sometimes M/D/Y,
        # so accomodate for this
        try:
            df['time'] = [dt.datetime.strptime(timestamp, '%Y-%m-%d').date() for timestamp in df['time']]
        except:
            df['time'] = [dt.datetime.strptime(timestamp, '%m/%d/%Y').date() for timestamp in df['time']]

        # extract unique dates and note number of time steps, and also create a date string
        # corresponding to these dates
        dates = np.sort(np.unique(df['time'].values))
        Nt = len(dates)
        datestrings = [dt.datetime.strftime(date,'%Y-%m-%d') for date in dates]

        ##########################################################################################
        # Intialize a master data frame to store all of the data by group (A, B, C, ...)
        ##########################################################################################
        
        # initialize a list that will contain all the column names for the master data frame, 
        # starting with the variables that are the same across all the parameters
        column_names = ['group','time','Volume (m^3)','Area (m^2)',
                                 param_name + ',' + 'Mass (Mg)', 
                                 param_name + ',' + 'dMass/dt (Mg/d)',
                                 param_name + ',' + 'Net Flux In (Mg/d)',
                                 param_name + ',' + 'Net Load (Mg/d)',
                                 param_name + ',' + 'Net Transport In (Mg/d)',
                                 param_name + ',' + 'Net Reaction (Mg/d)',
                                 param_name + ',' + 'dMass/dt, Balance Check (Mg/d)']

        # append fluxes from N, S, E, W, etc.
        for direction in direction_list:
            column_names.append(param_name + ',' + 'Flux In from %s (Mg/d)' % direction)
                    
        # append the reaction terms, which are unique to each parameter, also save the 
        # list of reaction names . this code assumes that the reaction terms are all the 
        # column names with a comma in them followed by a 'd'
        reaction_name_list = []
        for rt in df.columns:
            if ',d' in rt:
                column_names.append(rt + ' (Mg/d)')
                reaction_name_list.append(rt + ' (Mg/d)')
        
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
            
            # initialize total flux through this group/face pair at zero 
            flux_param = np.zeros(Nt, dtype=float)
            
            # loop through pairs of control volumes comprising total flux through this group/face, 
            # and add up the fluxes
            for flux_pair in flux_pairs_dict[flux_key]:  
                
                # each flux_pair is a 2-item list where the first item is the "from" 
                # control volume/polygon and the second item is the "to" control volume/polygon
                from_poly = flux_pair[0]
                to_poly = flux_pair[1]        
                
                # using logical indexing, extract the rows of the balance table data frames corresponding to 
                # the "from" polygon in this particular flux pair for each parameter, storing the values in a list
                rows_from_param = df.loc[df['Control Volume'] == ('polygon%d' % from_poly)]
                
                # only need one row to identify which control volume is the "to" polygon, so pick out the first row
                row_to_param = rows_from_param.iloc[0]
                    
                # count the number of fields labeled "To_polyN" where N is an integer
                n2poly = 0
                for key in row_to_param.keys():
                    if 'To_poly' in key:
                        n2poly = n2poly + 1
                
                # note that the way Zhenlin handled the "to" control volumes is by creating eight fields called 
                # "To_polyX" where X = 0, 1, ..., 7, that map to the numbers of all the "to" control volumes for each
                # "from" polygon. using our one example row, find the number X of the "to" polygon in the corresponding 
                # to the "to" control volume in our flux pair -- should be the same for all substances but doing them 
                # separately to narrow down origin of any potential errors
                polynum_param = None
                for p in range(0,n2poly):
                    if row_to_param['To_poly%d' % p] == to_poly:
                        polynum_param = p
                
                # if we can't find the "to" polygon, polynum will remain None, and the 
                # following will throw an assertion error. if you get this error, it means
                # there is either an error in the definition of the flux dictionary or an error
                # in Zhenlin's data tables, probably the dictionary
                assert(polynum_param+1)
                
                # add the flux for this control volume pair flux to the total flux for the group/side pair
                flux_param += rows_from_param["Flux%d" % polynum_param].values.copy()
                
            # add the total flux for this group/side pair to the master data frame in the appropriate place
            # switching direction because while Zhenlin's "To_poly" implies these fluxes are out, they are
            # actually in. divide by 1e6 to converg g/d to Mg/d
            ind_m = df_master['group'] == group
            df_master.loc[ind_m,param_name + ',Flux In from %s (Mg/d)' % NSEW] = flux_param.copy()/1.0e6
               
        ########################################################################################
        # For each group (A, B, C, ...) add up all the contributions from the multiple
        # control volumes that comprise the group to get the total area, volume, rate of 
        # change of mass, loads, and reaction rates
        ########################################################################################
        
        # loop through the groups
        for group in group_dict.keys():
    
            # initialize the area and volume of the group, 
            # which are the same for all the parameters, at zero
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
               
                # add volumes and areas
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
            
            # add volumes and areas
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
        # Add up the fluxes and reactions for this parameter to get the net and add to dataframe
        ########################################################################################
                    
        # add N,S,E,W,NE,etc components to get the net flux out and set value in the master data
        # frame
        df_master[param_name + ',Net Flux In (Mg/d)'] = 0
        for direction in direction_list:
            df_master[param_name + ',Net Flux In (Mg/d)'] += df_master[param_name + ',Flux In from %s (Mg/d)' % direction].values.copy()
        
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
        
        df_master.to_csv(os.path.join(output_dir,run_folder,'%s_Table_By_Group.csv' % param_name.lower()),index=False)
        
        
    



