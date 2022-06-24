
'''
Starting with the composite parameter balance tables generated in step2, group together some of the
reactions to create composite reactions. This is an important step where we identify groups of reacitons
that should cancel out to zero because they represent passing of mass from one substance to another (such
as nitrification passing N from NH4 to NO3, in the DIN budget), and thus we may identify mass conservation
errors. Allie April 2022
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

# run id 
#runid = 'FR16_28'
#runid = 'FR13_003'
#runid = 'FR13_023'
#runid = 'FR17_003'
#runid = 'FR17_017'
runid = 'FR18_005'

# is delta run?
is_delta = False

# path to base level balance tables (organized by runid within this folder)
balance_table_path = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables'

# call out any sets of reactions you want to sum together to appear as one reaction 
# in the composite parameter mass budgets, and give each composite reaction a name ...
composite_reaction_dict = {
     'Zoopl' : { 'Zoopl_VER_dZ_VERmor' : ['Zoopl_V,dZ_Vmor','Zoopl_E,dZ_Emor','Zoopl_R,dZ_Rmor'],
                 'Zoopl_VR_VRgr' : ['Zoopl_V,dZ_Vgr', 'Zoopl_R,dZ_Rgr']},
     'Algae' : {'Algae,dZ_Algae' : ['Diat,dZ_Diat', 'Green,dZ_Grn'],
                'Algae,dSedAlgae' : ['Diat,dSedDiat', 'Green,dSedGreen'], 
                'Algae,dPPAlgae' : ['Diat,dPPDiat','Green,dPPGreen'],
                'Algae,dcPPAlgae' : ['Diat,dcPPDiat','Green,dcPPGreen'],
                'Algae,dMrtAlgae' : ['Diat,dMrtDiat','Green,dMrtGreen'] },
     'DIN' : {
             'NH4,dMinDetN' : ['NH4,dMinDetNS1', 'NH4,dMinDetNS2'],
             'DIN,dDINUpt'  : ['NO3,dNO3Upt', 'NH4,dNH4Upt'],       
             'NO3,dDenit'   : ['NO3,dDenitWat', 'NO3,dDenitSed'],
             'ZERO: NO3,dNiDen'     : ['NO3,dNiDen'],                             # this is zero for the current model (because no sediment model??? not sure why)
             'ZERO: NH4,dNitrif + NO3,dNitrif'  : ['NO3,dNitrif', 'NH4,dNitrif'],     # these should cancel out    
             },
    'TN' : { 'NO3,dDenit' : ['NO3,dDenitWat', 'NO3,dDenitSed'], # this should be a SINK for TN
             'Algae,dSedAlgae' : ['Diat,dSedDiat', 'Green,dSedGreen'],  # this should be a SINK for TN
             'PON1,dSedPON1' : ['PON1,dSedPON1'], # this should be a SINK for TN 
             'NH4,dMinDetNS' : ['NH4,dMinDetNS1', 'NH4,dMinDetNS2'], # this should be a SOURCE for TN
             'NH4,dClam_NRes' : ['NH4,dM_NRes','NH4,dG4_NRes'], # this is a SOURCE for TN (clam pee)
             'PON1,dClam_NDef' : ['PON1,dM_NDef','PON1,dG4_NDef'], # this is a SOURCE for TN (clam poo)
             'Algae,dClam_Algae' : ['Diat,dM_Diat','Diat,dG4_Diat','Green,dM_Green','Green,dG4_Green'], # this is a SINK for TN (clams eat algae)
             'PON1,dClam_PON1' : ['PON1,dM_PON1','PON1,dG4_PON1'], # this is a SINK for TN (clams eat dead algae?)
             'ZERO: NO3,dNiDen'     : ['NO3,dNiDen'],                             # this is zero for the current model (because no sediment model??? not sure why)
             'ZERO: PON1,dCnvPPON1' : ['PON1,dCnvPPON1'],
             'ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN' : [
                                          'PON1,dResS1DetN', 
                                          'PON1,dResS2DetN', 
                                          'PON1,dResS1DiDN', 
                                          'PON1,dResS2DiDN'],   
             'ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet' : [
                                          'PON1,dZ_PON1',
                                          'PON1,dZ_NMrt',      
                                          'PON1,dM_NMrt',
                                          'PON1,dG4_NMrt',
                                          'PON1,dM_NSpDet',
                                          'PON1,dG4_NSpDet'],    
             'ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor' : [                                        
                                         'Zoopl_V,dZ_Vmor',  
                                          'Zoopl_R,dZ_Rmor',  
                                          'Zoopl_E,dZ_Emor'], 
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt' : 
                                                                                     ['Diat,dPPDiat',
                                                                                      'Diat,dcPPDiat',
                                                                                      'Green,dPPGreen',
                                                                                      'Green,dcPPGreen',
                                                                                      'NH4,dNH4Upt',
                                                                                      'NO3,dNO3Upt'],
             # identify groups of terms that sum to zero
             'ZERO: NO3,dNitrif + NH4,dNitrif' : ['NO3,dNitrif','NH4,dNitrif'],      # these two sum to zero
             'ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1' : ['DON,dCnvDPON1','PON1,dCnvDPON1'], # these two sum to zero
             'ZERO: NH4,dMinPON1 + PON1,dMinPON1' : ['NH4,dMinPON1','PON1,dMinPON1'],   # these two sum to zero
             'ZERO: NH4,dMinDON + DON,dMinDON' : ['NH4,dMinDON','DON,dMinDON'],       # these two sum to zero

             'ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut' :  [ 
                                        'Diat,dMrtDiat', 
                                        'Green,dMrtGreen',
                                        'PON1,dMortDetN',
                                        'NH4,dNH4Aut'], 
             'ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec' : [
                                         'NH4,dZ_NRes',
                                         'PON1,dZ_NDef',
                                         'PON1,dZ_NSpDet',
                                         'Diat,dZ_Diat',
                                         'Green,dZ_Grn',
                                         'Zoopl_V,dZ_Vgr',
                                         'Zoopl_R,dZ_SpwDet',
                                         'Zoopl_R,dZ_Rgr',
                                         'Zoopl_E,dZ_Ea',
                                         'Zoopl_E,dZ_Ec'
                                         ]
            },
    'TN_include_sediment' : {
             'NO3,dDenit' : ['NO3,dDenitWat','NO3,dDenitSed'], # this should be a SINK for TN
             'DetNS2,dBurS2DetN' : ['DetNS2,dBurS2DetN'], # this appears to act as a true SINK for TN as well
             'ZERO: NO3,dNiDen'     : ['NO3,dNiDen'],                             # this is zero for the current model (because no sediment model??? not sure why)
             'ZERO: PON1,dCnvPPON1' : ['PON1,dCnvPPON1'],
             'ZERO: PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN' : [
                                          'PON1,dResS1DetN', 
                                          'PON1,dResS2DetN', 
                                          'PON1,dResS1DiDN', 
                                          'PON1,dResS2DiDN'],   
             'ZERO: PON1,dZ_PON1 + PON1,dZ_NMrt + PON1,dM_NMrt + PON1,dG4_NMrt + PON1,dM_NSpDet + PON1,dG4_NSpDet' : [
                                          'PON1,dZ_PON1',
                                          'PON1,dZ_NMrt',      
                                          'PON1,dM_NMrt',
                                          'PON1,dG4_NMrt',
                                          'PON1,dM_NSpDet',
                                          'PON1,dG4_NSpDet'],    
             'ZERO: Zoopl_V,dZ_Vmor + Zoopl_R,dZ_Rmor' : [                                        
                                 'Zoopl_V,dZ_Vmor',  
                                  'Zoopl_R,dZ_Rmor',  
                                  'Zoopl_E,dZ_Emor'], 
             'ZERO: DetNS1,dZ_NMrtS1 + DetNS1,dZ_DNS1 + DetNS1,dM_DNS1 + DetNS1,dG4_DNS1' :   [
                         'DetNS1,dZ_NMrtS1',
                         'DetNS1,dZ_DNS1',
                         'DetNS1,dM_DNS1',
                         'DetNS1,dG4_DNS1'],
             'ZERO: DetNS1,dSWMinDNS1 + DetNS1,dResS1DetN + DetNS1,dSWBuS1DtN + DetNS1,dDigS1DetN' : [
                         'DetNS1,dSWMinDNS1',
                         'DetNS1,dResS1DetN',
                         'DetNS1,dSWBuS1DtN',
                         'DetNS1,dDigS1DetN'],
             'ZERO: DetNS2,dSWMinDNS2 + DetNS2,dResS2DetN + DetNS2,dDigS1DetN + DetNS2,dDigS2DetN' : [
                         'DetNS2,dSWMinDNS2',
                         'DetNS2,dResS2DetN',
                         'DetNS2,dDigS1DetN',
                         'DetNS2,dDigS2DetN'],
             'ZERO: Mussel_R,dM_SpwDet + Grazer4_R,dG4_SpwDet' : ['Mussel_R,dM_SpwDet','Grazer4_R,dG4_SpwDet'],
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt' : 
                                                                                     ['Diat,dPPDiat',
                                                                                      'Diat,dcPPDiat',
                                                                                      'Green,dPPGreen',
                                                                                      'Green,dcPPGreen',
                                                                                      'NH4,dNH4Upt',
                                                                                      'NO3,dNO3Upt'],
             # identify groups of terms that sum to zero
             'ZERO: NO3,dNitrif + NH4,dNitrif' : ['NO3,dNitrif','NH4,dNitrif'],      # these two sum to zero
             'ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1' : ['DON,dCnvDPON1','PON1,dCnvDPON1'], # these two sum to zero
             'ZERO: NH4,dMinPON1 + PON1,dMinPON1' : ['NH4,dMinPON1','PON1,dMinPON1'],   # these two sum to zero
             'ZERO: NH4,dMinDON + DON,dMinDON' : ['NH4,dMinDON','DON,dMinDON'],       # these two sum to zero
             'ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut' :  [ 
                                        'Diat,dMrtDiat', 
                                        'Green,dMrtGreen',
                                        'PON1,dMortDetN',
                                        'NH4,dNH4Aut'], 
             'ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec' : [
                                 'NH4,dZ_NRes',
                                 'PON1,dZ_NDef',
                                 'PON1,dZ_NSpDet',
                                 'Diat,dZ_Diat',
                                 'Green,dZ_Grn',
                                 'Zoopl_V,dZ_Vgr',
                                 'Zoopl_R,dZ_SpwDet',
                                 'Zoopl_R,dZ_Rgr',
                                 'Zoopl_E,dZ_Ea',
                                 'Zoopl_E,dZ_Ec'],
             'ZERO: PON1,dSedPON1 + DetNS1,dSedPON1' : ['PON1,dSedPON1',
                                                        'DetNS1,dSedPON1'],
             'ZERO: NH4,dMinDetNS1 + NH4,dMinDetNS2 + DetNS1,dMinDetNS1 + DetNS2,dMinDetNS2' : ['NH4,dMinDetNS1', 
                                                                                                'NH4,dMinDetNS2', 
                                                                                                'DetNS1,dMinDetNS1',
                                                                                                'DetNS2,dMinDetNS2'],
             'ZERO: Mussel_V,dM_Vmor + Mussel_E,dM_Emor + Mussel_R,dM_Rmor + DetNS1,dM_NMrtS1' : ['Mussel_V,dM_Vmor',
                                                                                                     'Mussel_E,dM_Emor',
                                                                                                     'Mussel_R,dM_Rmor', 
                                                                                                     'DetNS1,dM_NMrtS1'],
             'ZERO: Grazer4_V,dG4_Vmor + Grazer4_E,dG4_Emor + Grazer4_R,dG4_Rmor + DetNS1,dG4_NMrtS1' : ['Grazer4_V,dG4_Vmor',
                                                                                                            'Grazer4_E,dG4_Emor',
                                                                                                            'Grazer4_R,dG4_Rmor',
                                                                                                            'DetNS1,dG4_NMrtS1'],
             'ZERO: NH4,dM_NRes + PON1,dM_NDef + PON1,dM_PON1 + Diat,dM_Diat + Mussel_V,dM_Vgr + Mussel_E,dM_Ea + Mussel_E,dM_Ec + Mussel_R,dM_Rgr' : ['NH4,dM_NRes',
                       'PON1,dM_NDef',
                       'PON1,dM_PON1',
                       'Diat,dM_Diat',
                       'Mussel_V,dM_Vgr',
                       'Mussel_E,dM_Ea',
                       'Mussel_E,dM_Ec',
                       'Mussel_R,dM_Rgr'],
             'ZERO: NH4,dG4_NRes + PON1,dG4_NDef + PON1,dG4_PON1 + Diat,dG4_Diat + Grazer4_V,dG4_Vgr + Grazer4_E,dG4_Ea + Grazer4_E,dG4_Ec + Grazer4_R,dG4_Rgr' : ['NH4,dG4_NRes',
                         'PON1,dG4_NDef',
                         'PON1,dG4_PON1',
                         'Diat,dG4_Diat',
                         'Grazer4_V,dG4_Vgr',
                         'Grazer4_E,dG4_Ea',
                         'Grazer4_E,dG4_Ec',
                         'Grazer4_R,dG4_Rgr'],
             'ZERO: DetNS1,dBurS1DetN + DetNS2,dBurS1DetN' : ['DetNS1,dBurS1DetN', 'DetNS2,dBurS1DetN'],
             'ZERO: Diat,dSedDiat + Green,dSedGreen + DetNS1,dSedAlgN' :  ['Diat,dSedDiat','Green,dSedGreen','DetNS1,dSedAlgN'] 
    }
}

# float format for csv files
float_format = '%1.16e'

# if is delta, add "Delta_" prefix to runid
if is_delta:
    runid = 'Delta_' + runid

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

# define a function that returns a list of columns that are not reactions
def get_not_rx_list(df, rx_list):
    not_rx_list = []
    columns = df.columns
    for column in columns:
        if not column in rx_list:
            not_rx_list.append(column)
    return not_rx_list


###################################################################
# MAIN
###################################################################

# run folder is run id within balance table folder
run_folder = os.path.join(balance_table_path,runid)

# setup logging to file and to screen 
try:
    logging.basicConfig(
        level=logging.INFO,
        mode='w',
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(run_folder,"log_step3.log")),
            logging.StreamHandler(sys.stdout)
        ]
    )
except:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(run_folder,"log_step3.log")),
            logging.StreamHandler(sys.stdout)
        ]
    )

# add some basic info to log file
user = os.getlogin()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('Composite balance tables with "Grouped Rx" were produced on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))

# loop through the composite parameters
for composite_param in composite_reaction_dict.keys():

    # log
    logging.info('Grouping reaction terms for %s' % composite_param)

    # input and output table names and paths
    input_table_name = composite_param.lower() + '_Table_Ungrouped_Rx.csv'
    input_table_path = os.path.join(run_folder, input_table_name)
    output_table_name = composite_param.lower() + '_Table.csv'
    output_table_path = os.path.join(run_folder, output_table_name)
    
    # get dictionary of compositie reactions
    composite_reactions = composite_reaction_dict[composite_param]

    #  read input table, which is the balance table for the composite parameters including every single reaction, totally ungrouped
    df_composite = pd.read_csv(input_table_path)

    # get the list of columns that contain reactions, and the list that don't
    rx_list = get_rx_list(df_composite, composite_param)
    not_rx_list = get_not_rx_list(df_composite, rx_list)

    # do some logging of progress
    logging.info('    The %s budget includes the following reaction terms:' % composite_param)
    for rx in rx_list:
        logging.info('        %s' % rx)
    logging.info('    Now we will group some of these reactions together, taking care not to lose track of any of them ...')

    # initialize the balance table with the grouped parameters
    df_grouped = df_composite[not_rx_list].copy(deep=True)

    # now add all the grouped reactions to the balance table, removing their components from the list of 
    # all the component reactions as we go along
    grouped_rx_list = list(composite_reactions.keys())
    remaining_rx_list = rx_list.copy()
    for grouped_rx in grouped_rx_list:

        # log
        logging.info('        Summing up componenets of %s...' % grouped_rx)
        
        # initialize the grouped reaction at zero
        df_grouped[grouped_rx] = 0

        # sum up its components
        component_reactions = composite_reactions[grouped_rx]
        for rx in component_reactions:

            # log 
            logging.info('        ... %s' % rx)
            
            # if the reaction isn't in the original list of reactions, skip it
            if not (rx in rx_list):
                logging.info('            %s not found, probably not included in this model run, skipping it...' % rx)
                continue

            # if the reaction WAS in the original list, but it isn't in the list of remaining reactions,
            # it must have been double counted, so raise an exception
            elif not rx in remaining_rx_list:
                raise Exception('%s is double counted in composite reactions for %s' % (rx, composite_param))
            
            # add the contribution of the component reaction to its grouped reaction
            df_grouped[grouped_rx] = df_grouped[grouped_rx].values + df_composite[rx].values # add contribution to grouped reaction
            
            # remove this reaction from the list of remaining reactions
            remaining_rx_list.remove(rx)

    # now be sure to include all the remaining reactions
    logging.info('        The following reactions were not used in any grouped reactions, adding them to the final budget...')
    for rx in remaining_rx_list:
        logging.info('        ... %s' % rx)
        df_grouped[rx] = df_composite[rx].values

    # double check the final reaction list
    rx_list = get_rx_list(df_grouped, composite_param) 
    logging.info('    The %s budget includes the following reaction terms after grouping:' % composite_param)
    for rx in rx_list:
        logging.info('        %s' % rx)

    # now save the composite balance table
    logging.info('    Saving composite parameter table with grouped reactions %s' % output_table_path)
    df_grouped.to_csv(output_table_path, index=False, float_format=float_format)

