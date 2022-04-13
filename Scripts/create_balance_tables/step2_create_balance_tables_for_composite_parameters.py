
'''
Starting with raw parameters found in 'Balance_Table_All_Parameters_By_Group.csv' file, 
calculates fluxes and rates for composite parameters defined in the user input section by
constituents and reaction groupings.
'''

#################################################################
# IMPORT PACKAGES
#################################################################

import sys, os, re
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import logging
import datetime
import socket
hostname = socket.gethostname()
if hostname == 'richmond':
    raise Exception("Do not run this script on richmond until we update the conda environment... run on chicago or your laptop instead")


#################################################################
# USER INPUT
#################################################################

# open bay, has clams: G141_13to18_132

# run id 
runid = 'FR13_003'
#runid = 'Delta_FR16_28'

# path to the lsp file
#lsp_path = r'W:\open_bay\BGC_model\Full_res\WY2013\FR13_021\sfbay_dynamo000.lsp'
#lsp_path = r'X:\hpcshared\Delta\BGC_model\Full_res\FR16_28\FR16_28\wy2016.lsp'
#lsp_path = '/richmondvol1/hpcshared/Delta/BGC_model/Full_res/FR16_28/FR16_28/wy2016.lsp'
lsp_path = '/richmondvol1/hpcshared/Full_res/WY2013/FR13_003/sfbay_dynamo000.lsp'


# path to base level balance tables (organized by runid within this folder)
#balance_table_path = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Balance_Tables'
balance_table_path = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables'

# list of composite parameters -- include all possible component parameters, script will automatically 
# check if they were included in this run and leave it out if needed
composite_parameters = {
    'DIN' : ['NH4', 'NO3'],
    'TN'  : ['NH4', 'NO3', 'PON1', 'DON', 'Diat', 'Green', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R'], 
    'TP'  : ['PO4', 'POP1', 'DOP', 'Diat', 'Green', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'DetNS' : ['DetNS1', 'DetNS2'],
    'DetPS' : ['DetPS1', 'DetPS2'],
    'DetSi' : ['DetSiS1', 'DetSiS2'],
    'Zoopl' : ['Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'Grazer4' : ['Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'Mussel' : ['Zoopl_V', 'Zoopl_E', 'Zoopl_R'], 
    'TP_include_sediment' : ['PO4', 'POP1', 'DOP', 'DetPS1', 'DetPS2', 
                             'Diat', 'Green', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R',
                             'Mussel_V','Mussel_E','Mussel_R','Grazer4_V','Grazer4_E','Grazer4_R'],
    'TN_include_sediment' : ['NH4', 'NO3', 'PON1', 'DON', 'DetNS1', 'DetNS2', 
                             'Diat', 'Green', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R',
                             'Mussel_V','Mussel_E','Mussel_R','Grazer4_V','Grazer4_E','Grazer4_R']
}

# is the budget of each composite parameter in grams of C, N, P, etc?
composite_bases = {
    'DIN' : 'N',
    'TN'  : 'N', 
    'TP'  : 'P', 
    'TN_include_sediment'  : 'N', 
    'TP_include_sediment'  : 'P', 
    'DetNS' : 'N',
    'DetPS' : 'P',
    'DetSi' : 'Si',
    'Zoopl' : 'C',
    'Grazer4' : 'C',
    'Mussel' : 'C'
}

# float format for csv files
float_format = '%1.16e'

#################################################################
# FUNCTIONS
#################################################################

# define a function that returns a list of the reactions included in a balance table
def get_rx_list(df, substance):    
    rx_list = []
    columns = df.columns  
    for column in columns:
        if not column in ['time', 'Control Volume', 'Concentration (mg/l)', 'Volume',
                          'Volume (Mean)', 'Area', 'dVar/dt', '%s,Loads in' % substance, '%s,Loads out' % substance,
                          '%s,Transp in' % substance, '%s,Transp out' % substance]:
            if not 'Flux' in column:
                if not 'To_poly' in column:
                    rx_list.append(column)
    return rx_list

###################################################################
# MAIN
###################################################################

# run folder is run id within balance table folder
run_folder = os.path.join(balance_table_path,runid)

# setup logging to file and to screen 
logging.basicConfig(
    level=logging.INFO,
    mode='w',
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(run_folder,"log_step2.log")),
        logging.StreamHandler(sys.stdout)
    ]
)

# add some basic info to log file
user = os.getlogin()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('Composite balance tables with "Ungrouped Rx" were produced on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))

# read the lsp file and return a list of the substances
substances = []
with open(lsp_path,'r') as f:
    for line in f.readlines():
        if '-fluxes for [' in line:
            # the following ugly expression returns the string between the brackets, stripped of whitespace:
            substances.append(line[line.find("[")+1:line.find("]")].strip()) 
logging.info('The following substances were found in %s:' % lsp_path)
for substance in substances:
    logging.info('    %s' % substance)

# scan the lsp file for mentions of N:C ratios and P:C ratios
NC_ratios = {}
PC_ratios = {}
with open(lsp_path,'r') as f:
    lines = f.readlines()
    for i in range(len(lines)-1):
        line1 = lines[i]
        line2 = lines[i+1]
        if 'N:C ratio Diatoms' in line1:
            NC_ratios['Diat'] = float(line2.split(':')[1])
        elif 'N:C ratio Greens' in line1:
            NC_ratios['Green'] = float(line2.split(':')[1])
        elif 'N:C ratio of DEB Zooplankton' in line1:
            NC_ratios['Zoopl_V'] = float(line2.split(':')[1])
            NC_ratios['Zoopl_E'] = float(line2.split(':')[1])
            NC_ratios['Zoopl_R'] = float(line2.split(':')[1])
        elif 'N:C ratio of DEB Grazer4' in line1:
            NC_ratios['Grazer4_V'] = float(line2.split(':')[1])
            NC_ratios['Grazer4_E'] = float(line2.split(':')[1])
            NC_ratios['Grazer4_R'] = float(line2.split(':')[1])
        elif 'N:C ratio of DEB Mussel' in line1:
            NC_ratios['Mussel_V'] = float(line2.split(':')[1])
            NC_ratios['Mussel_E'] = float(line2.split(':')[1])
            NC_ratios['Mussel_R'] = float(line2.split(':')[1])
        if 'P:C ratio Diatoms' in line1:
            PC_ratios['Diatom'] = float(line2.split(':')[1])
        elif 'P:C ratio Greens' in line1:
            PC_ratios['Green'] = float(line2.split(':')[1])
        elif 'P:C ratio of DEB Zooplankton' in line1:
            PC_ratios['Zoopl_V'] = float(line2.split(':')[1])
            PC_ratios['Zoopl_E'] = float(line2.split(':')[1])
            PC_ratios['Zoopl_R'] = float(line2.split(':')[1])
        elif 'P:C ratio of DEB Grazer4' in line1:
            PC_ratios['Grazer4_V'] = float(line2.split(':')[1])
            PC_ratios['Grazer4_E'] = float(line2.split(':')[1])
            PC_ratios['Grazer4_R'] = float(line2.split(':')[1])
        elif 'P:C ratio of DEB Mussel' in line1:
            PC_ratios['Mussel_V'] = float(line2.split(':')[1])
            PC_ratios['Mussel_E'] = float(line2.split(':')[1])
            PC_ratios['Mussel_R'] = float(line2.split(':')[1])
logging.info('The following N:C ratios were found in %s:' % lsp_path)
keys = list(NC_ratios.keys())
keys.sort()
for key in keys:
    logging.info('    %s : %f' % (key, NC_ratios[key]))
logging.info('The following P:C ratios were found in %s:' % lsp_path)
keys = list(PC_ratios.keys())
keys.sort()
for key in keys:
    logging.info('    %s : %f' % (key, PC_ratios[key]))

# load all the balance tables into a dictionary
logging.info('Reading balance tables from %s:' % run_folder)
df_dict = {}
for substance in substances:
    table_name = '%s_Table.csv' % substance.lower()
    logging.info('    %s' % table_name)
    df = pd.read_csv(os.path.join(run_folder,table_name))
    df_dict[substance] = df.copy(deep=True)    

# find the number of adjacent polygons included in the fluxes, using the NH4 balance table
logging.info('Using the NH4 balance table to check number of adjacent polygons included in the fluxes...')
polymax = 0
for column in df_dict['NH4'].columns:
    if 'To_poly' in column:
        n = int(column.strip('To_poly'))
        if n>polymax:
            polymax = n
logging.info('Tracking %d fluxes to adacent polygons' % (polymax+1))

# make a base balance table from the NH4 balance tables, to reuse when initializing composite parameter tables
# ... these columns have the same values in all balance tables
df_base = df_dict['NH4'][['time', 'Control Volume', 'Volume',
       'Volume (Mean)', 'Area']].copy(deep=True)
# ... initialize these columns at zero
df_base['Concentration (mg/l)'] = 0.0
df_base['dVar/dt'] = 0.0
for n in range(polymax+1):
    df_base['To_poly%d' % n] = df_dict['NH4']['To_poly%d' % n].values
    df_base['Flux%d' % n] = 0.0

# make a list of the flux columns to track when adding up the component parameters
flux_list = []
for n in range(polymax+1):
    flux_list.append('Flux%d' % n)

# make a list of the loading and transport terms, with an %s placeholder for the substance
load_transp_list = ['%s,Loads in','%s,Loads out','%s,Transp in','%s,Transp out']

# loop through the composite parameters
for composite_param in composite_parameters.keys():
    
    # log
    logging.info('Deriving balance table for %s' % composite_param)

    # check the base element for the composite parameter budget
    composite_base = composite_bases[composite_param]
    logging.info('    Base element is %s' % composite_base)

    # initialize the balance table
    df_composite = df_base.copy(deep=True)

    # initialize columns for the loading and transport
    for column in load_transp_list:
        df_composite[column % composite_param] = 0

    # loop through the component parameters making up the composite parameter
    for component_param in composite_parameters[composite_param]:

        # log component parameter
        logging.info('    Adding %s budget ...' % component_param)

        # check if this parameter was included in the model run, and if not, skip it
        if not component_param in substances:
            logging.info('       %s not found in list of substances, skipping this component of %s' % (component_param, composite_param))
            continue    

        # figure out the correct stoichiometric multiplier
        multiplier = 1
        if composite_base == 'N':
            if component_param in NC_ratios.keys():
                multiplier = NC_ratios[component_param]
        elif composite_base == 'P':
            if component_param in PC_ratios.keys():
                multiplier = PC_ratios[component_param]
        logging.info('        Using the following stoichiometric multiplier for %s: %f' % (component_param, multiplier))

        # add up the concentrations, the change in mass, and the fluxes
        for column in (['Concentration (mg/l)', 'dVar/dt'] + flux_list):
            logging.info('        Adding contribution to %s' % column)
            df_composite[column] = df_composite[column].values + multiplier * df_dict[component_param][column].values

        # add up the loadings and transport terms
        for column in load_transp_list:
            logging.info('        Adding contribution to %s' % (column % composite_param))
            df_composite[column % composite_param] = df_composite[column % composite_param].values + multiplier * df_dict[component_param][column % component_param].values

        # add the reactions
        rx_list = get_rx_list(df_dict[component_param], component_param)
        for column in rx_list:
            logging.info('        Adding the reaction term %s' % column)
            df_composite[column] = multiplier * df_dict[component_param][column].values

    # now save the composite balance table
    table_name = composite_param.lower() + '_Table_Ungrouped_Rx.csv'
    table_path = os.path.join(run_folder, table_name)
    logging.info('    Saving composite parameter table %s' % table_path)
    df_composite.to_csv(table_path, index=False, float_format=float_format)

