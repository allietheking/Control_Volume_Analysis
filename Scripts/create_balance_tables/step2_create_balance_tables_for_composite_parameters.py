
'''
Starting with balance tables for modeled parameters, create balance tables for composite parameters
such as DIN and TN. Automatically finds correct stoichiometric mutipliers from *.lsp file. All the 
reactions in each of the component parameters are included separately. They may be consolidated into 
composite reaction terms in "step3". Allie April 2022
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

# if running the script alone, load the configuration module (in this folder)
if __name__ == "__main__":

    import importlib
    import step0_config
    importlib.reload(step0_config)

#################################################################
# USER INPUT
#################################################################

# get variables out of the configuration module (see step0_config.py in this folder)
from step0_config import runid, is_delta, balance_table_dir, model_inout_dir, float_format

# ugly complicated process to get the water year from the runid
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

# path to the lsp file
#lsp_path = '/richmondvol1/hpcshared/Delta/BGC_model/Full_res/FR16_28/FR16_28/wy2016.lsp'
#lsp_path = '/hpcvol2/open_bay/BGC_model/Full_res/%s/%s/sfbay_dynamo000.lsp' % (water_year,runid)
if 'FR' in runid:
    lsp_path = os.path.join(model_inout_dir,'Full_res',water_year,runid,'sfbay_dynamo000.lsp')
else:
    lsp_path = os.path.join(model_inout_dir,'Grid141',water_year,runid,'sfbay_dynamo000.lsp')

# list of composite parameters -- include all possible component parameters, script will automatically 
# check if they were included in this run and leave it out if needed
composite_parameters = {
    'N-Algae' : ['Diat', 'Green','DiatS1'],
    'N-Zoopl' : ['Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'Algae' : ['Diat', 'Green','DiatS1'],
    'Zoopl' : ['Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'DIN' : ['NH4', 'NO3'],
    'TN'  : ['NH4', 'NO3', 'PON1', 'DON', 'Diat', 'DiatS1', 'Green', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R'], 
    'TN_include_sediment' : ['NH4', 'NO3', 'PON1', 'DON', 'DetNS1', 'DetNS2', 
                             'Diat', 'Green', 'DiatS1', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R',
                             'Mussel_V','Mussel_E','Mussel_R','Grazer4_V','Grazer4_E','Grazer4_R'],    
    'DetNS' : ['DetNS1', 'DetNS2'],                  
    'TP'  : ['PO4', 'POP1', 'DOP', 'Diat', 'Green', 'DiatS1', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
    'TP_include_sediment' : ['PO4', 'POP1', 'DOP', 'DetPS1', 'DetPS2', 
                             'Diat', 'Green', 'DiatS1', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R',
                             'Mussel_V','Mussel_E','Mussel_R','Grazer4_V','Grazer4_E','Grazer4_R'],
    'DetPS' : ['DetPS1', 'DetPS2'],
    'DetSi' : ['DetSiS1', 'DetSiS2'],
    'Grazer4' : ['Grazer4_V', 'Grazer4_E', 'Grazer4_R'],
    'Mussel' : ['Mussel_V', 'Mussel_E', 'Mussel_R'], 
    'Clams' : ['Grazer4_V', 'Grazer4_E', 'Grazer4_R', 'Mussel_V', 'Mussel_E', 'Mussel_R']
}

# is the budget of each composite parameter in grams of C, N, P, etc?
composite_bases = {
    'N-Algae' : 'N',
    'N-Zoopl' : 'N',
    'Algae' : 'C',
    'Zoopl' : 'C',
    'DIN' : 'N',
    'TN'  : 'N', 
    'TN_include_sediment'  : 'N', 
    'DetNS' : 'N',
    'TP'  : 'P', 
    'TP_include_sediment'  : 'P', 
    'DetPS' : 'P',
    'DetSi' : 'Si',
    'Grazer4' : 'C',
    'Mussel' : 'C',
    'Clams' : 'C'
}

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

def logger_cleanup():

    ''' check for open log files, close them, release handlers'''

    # clean up logging
    logger = logging.getLogger()
    handlers = list(logger.handlers)
    if len(handlers)>0:
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler) 

###################################################################
# MAIN
###################################################################

# setup logging to file and to screen 
logger_cleanup()
logging.basicConfig(
level=logging.INFO,
format="%(asctime)s [%(levelname)s] %(message)s",
handlers=[
    logging.FileHandler(os.path.join(balance_table_dir,"log_step2.log"),'w'),
    logging.StreamHandler(sys.stdout)
])

# add some basic info to log file
user = os.getlogin()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('Composite balance tables with "Ungrouped Rx" were produced on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))

# log configuration variables
logging.info('The following global variables were loaded from step0_conf.py:')
logging.info('    runid = %s' % runid)
logging.info('    is_delta = %r' % is_delta)
logging.info('    balance_table_dir = %s' % balance_table_dir)
logging.info('    model_inout_dir = %s' % model_inout_dir)
logging.info('    float_format = %s' % float_format)

# do some checks to make sure we're reading the right lsp file, makes some assumptiions about how files
# are organized on our servers, so if that changes may have to debug this part
if not runid in lsp_path:
    raise Exception('Danger: *.lsp path looks incorrect, aborting')
if is_delta:
    if (not (('Delta' in lsp_path) or ('delta' in lsp_path))):
        raise Exception('Danger: *.lsp path looks incorrect, aborting')

# scan the lsp file for mentions of N:C ratios and P:C ratios
NC_ratios = {}
PC_ratios = {}
with open(lsp_path,'r') as f:
    lines = f.readlines()
    for i in range(len(lines)-1):
        line1 = lines[i]
        line2 = lines[i+1]
        if 'NCRatDiat ' in line1:
            NC_ratios['Diat'] = float(line2.split(':')[1])
        elif 'N:C ratio Greens' in line1:
            NC_ratios['Green'] = float(line2.split(':')[1])
        elif 'NCRatDiatS ' in line1:
            NC_ratios['DiatS1'] = float(line2.split(':')[1])
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
        if 'PCRatDiat ' in line1:
            PC_ratios['Diat'] = float(line2.split(':')[1])
        elif 'P:C ratio Greens' in line1:
            PC_ratios['Green'] = float(line2.split(':')[1])
        elif 'PCRatDiatS ' in line1:
            PC_ratios['DiatS1'] = float(line2.split(':')[1])
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

# load all the balance tables into a dictionary
logging.info('Reading balance tables from %s:' % balance_table_dir)
df_dict = {}
for substance in substances.copy():
    table_name = '%s_Table.csv' % substance.lower()
    try:
        df = pd.read_csv(os.path.join(balance_table_dir,table_name))
    except:
        logging.info('    Did not find %s' % table_name)
        substances.remove(substance)
    else:
        logging.info('    Succesfully read %s' % table_name)
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

    # boolean to detect if any of the component parameters were included in the model run, if not, skip this composite parameter
    any_components = False

    # loop through the component parameters making up the composite parameter
    for component_param in composite_parameters[composite_param]:

        # log component parameter
        logging.info('    Adding %s budget ...' % component_param)

        # check if this parameter was included in the model run, and if not, skip it
        if not component_param in substances:
            logging.info('       %s not found in list of substances, skipping this component of %s' % (component_param, composite_param))
            continue    
        elif not component_param in df_dict.keys():
            logging.info('       ERROR: substance %s was modeled and is needed to compute %s, but could not find the %s_Table.csv file' % (component_param, composite_param, component_param))
            logging.info('       What to do? Add %s to substance_list in configuration file and re-run step1_create_balance_tables.py, then try again' % component_param)
            raise Exception('substance %s was modeled and is needed to compute %s, but could not find the %s_Table.csv file\n' % (component_param, composite_param, component_param) + 
                            'Add %s to substance_list in configuration file and re-run step1_create_balance_tables.py, then try again' % component_param)
        else:
            any_components = True

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
    table_path = os.path.join(balance_table_dir, table_name)
    if any_components:
        logging.info('    Saving composite parameter table %s' % table_path)
        df_composite.to_csv(table_path, index=False, float_format=float_format)
    else:
        logging.info('    No balance tables for component parameters of %s were found, so no composite table was created' % composite_param)

# clean up logging
logger_cleanup()