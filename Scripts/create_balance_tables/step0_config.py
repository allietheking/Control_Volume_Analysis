
''' 
Configuration module for ALL of the create_balance_table scripts 
This is where we set global variables that are used by all the scripts. 

alliek august 2022
'''

##############################
# user input
##############################

# this is the run you want to process
runid = 'G141_13to18_207'

# base level substances to process. set to string 'all' or a list of substance strings -- warning, processing 
# all of them takes a long time and uses a lot of space (this is used in step1_create_balance_tables.py)
# substance_list = 'all'
substance_list = ['continuity', 'nh4', 'no3', 'pon1', 'don', 'diat', 'diats1', 'green', 'oxy', 'zoopl_e', 'zoopl_r', 'zoopl_v', 'detns1', 'detns2', 'oons1', 'oons2']

# list of substances we think we are actually going to want to plot -- to save space, the 
# step5_compile_balance_tables_into_groups.py and step6_aggregate_in_time.py scripts will only process these substances
plot_substance_list = ['tn_include_sediment','tn','din','nh4','no3','don',
                        'pon1','n-algae','n-zoopl','totaldetns','detns1','detns2','oons1','oons2',
                        'oxy','algae','diat','green','diats1','zoopl']



# list of time averaging schemes to apply (saves space to skip some if we don't need them) 
# (this is used in step6_aggregate_in_time.py)
#tavg_list = ['Cumulative', 'Filtered', 'Annual', 'Seasonal', 'Monthly', 'Weekly']
tavg_list = ['Cumulative', 'Filtered', 'Seasonal']

# float format for csv files
float_format = '%1.6e'

# set tolerance for mass conservation error as a percentage of whatever variable we chose to normalize by (this is used in step4_check_mass_conservation.py)
error_tol_percent = 0.1

# abort if error tolerance is exceeded? set to False for troubleshooting only. (this is used in step4_check_mass_conservation.py)
abort_for_mass_cons_error = True

# delete all balance tables before re-running step1_create_balance_tables?
delete_balance_tables = True

# is this a delta run? a couple of notes about using these scripts for delta runs
# 1. for now we assume it is a full resolution delta run
# 2. since aggregated groups are not defined for the delta, do not try to run step5_compile_balance_tables_into_groups.py, 
#    but it is ok to run step6_aggregate_in_time.py without running step5_compile_balance_tables_into_groups.py
# 3. originally, i added the capacity to process delta runs to make sure mass of N is conserved, and for the most part it is, 
#    (there's a small leak when N is passed from DIN to algae because we're using an old version of DWAQ)
is_delta = False

# base directory for model input, namely the shapefiles (this definitely runs on linux, in theory can also run this in windows and use mounted drive)
#model_inout_dir = 'X:\hpcshared'
model_inout_dir = '/richmondvol1/hpcshared'

# stompy directory (stompy is called in step0_create_balance_tables.py to create dwaq_hist_bal.nc from the *-his.bal file if needed)
stompy_dir = '/opt/software/rusty/stompy/newest_commit/stompy'

# if the user sets the following variables to None, they are calculated automatically, by making some assumptions about
# how our computers are organized (see below)
poly_path = None
tran_path = None
group_def_path = None
group_con_path = None
run_dir = None
lsp_path = None
balance_table_dir = None

# for some parameters, the dwaq_hist.nc file does not have the correct units ... use this to override the units
# (this is used only in step1_create_balance_tables.py, where we compute dMass/dt from concentraitons in the *.his
# file, so we need to know if the concentraiton units are g/m3 or g/m3)
units_override = {'zoopl_e' : 'gC/m3',
                  'zoopl_r' : 'gC/m3',
                  'zoopl_v' : 'gC/m3',
                  'zoopl_n' : '#/m3',}

# list of composite parameters -- include all possible component parameters, script will automatically 
# check if they were included in this run and leave it out if needed (this is used in 
# step2_create_balance_tables_for_composite_parameters.py)
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
    'TotalDetNS' : ['DetNS1','OONS1','DetNS2','OONS2'],                  
#    'TP'  : ['PO4', 'POP1', 'DOP', 'Diat', 'Green', 'DiatS1', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R'],
#    'TP_include_sediment' : ['PO4', 'POP1', 'DOP', 'DetPS1', 'DetPS2', 
#                             'Diat', 'Green', 'DiatS1', 'Zoopl_V', 'Zoopl_E', 'Zoopl_R',
#                             'Mussel_V','Mussel_E','Mussel_R','Grazer4_V','Grazer4_E','Grazer4_R'],
#    'TotalDetPS' : ['DetPS1', 'DetPS2','OOPS1', 'OOPS2'],
#    'TotalDetSi' : ['DetSiS1', 'DetSiS2','OOSiS1', 'OOSiS2'],
#    'Grazer4' : ['Grazer4_V', 'Grazer4_E', 'Grazer4_R'],
#    'Mussel' : ['Mussel_V', 'Mussel_E', 'Mussel_R'], 
#    'Clams' : ['Grazer4_V', 'Grazer4_E', 'Grazer4_R', 'Mussel_V', 'Mussel_E', 'Mussel_R']
}

# is the budget of each composite parameter in grams of C, N, P, etc? (this is used in 
# step2_create_balance_tables_for_composite_parameters.py)
composite_bases = {
    'N-Algae' : 'N',
    'N-Zoopl' : 'N',
    'Algae' : 'C',
    'Zoopl' : 'C',
    'DIN' : 'N',
    'TN'  : 'N', 
    'TN_include_sediment'  : 'N', 
    'TotalDetNS' : 'N',
    'TP'  : 'P', 
    'TP_include_sediment'  : 'P', 
    'TotalDetPS' : 'P',
    'TotalDetSi' : 'Si',
    'Grazer4' : 'C',
    'Mussel' : 'C',
    'Clams' : 'C'
}

# call out any sets of reactions you want to sum together to appear as one reaction 
# in the composite parameter mass budgets, and give each composite reaction a name ...
# (this is used in step3_group_reactions_for_composite_parameters.py)
composite_reaction_dict = {
     'N-Algae' : {'Algae,dZ_Algae' : ['Diat,dZ_Diat', 'Green,dZ_Grn'],
                'Algae,dSedAlgae' : ['Diat,dSedDiat', 'Green,dSedGreen'], 
                'Algae,dPPAlgae' : ['Diat,dPPDiat','Green,dPPGreen','DiatS1,dPPDiatS1'],
                'Algae,dcPPAlgae' : ['Diat,dcPPDiat','Green,dcPPGreen'],
                'Algae,dMrtAlgae' : ['Diat,dMrtDiat','Green,dMrtGreen','DiatS1,dMrtDiatS1'],
                'DiatS1,dBurS1Diat' : ['DiatS1,dBurS1Diat'], # burial of benthic algae is currently a true sink in the model
                'ZERO: DiatS1,dResS1Diat' : ['DiatS1,dResS1Diat'], 
                'ZERO: DiatS1,dSWBuS1Dia' : ['DiatS1,dSWBuS1Dia'],
                'ZERO: DiatS1,dDigS1Diat' : ['DiatS1,dDigS1Diat'] },
     'N-Zoopl' : { 'Zoopl_VER_dZ_VERmor' : ['Zoopl_V,dZ_Vmor','Zoopl_E,dZ_Emor','Zoopl_R,dZ_Rmor'],
                 'Zoopl_VR_VRgr' : ['Zoopl_V,dZ_Vgr', 'Zoopl_R,dZ_Rgr']},
     'Algae' : {'Algae,dZ_Algae' : ['Diat,dZ_Diat', 'Green,dZ_Grn'],
                'Algae,dSedAlgae' : ['Diat,dSedDiat', 'Green,dSedGreen'], 
                'Algae,dPPAlgae' : ['Diat,dPPDiat','Green,dPPGreen','DiatS1,dPPDiatS1'],
                'Algae,dcPPAlgae' : ['Diat,dcPPDiat','Green,dcPPGreen'],
                'Algae,dMrtAlgae' : ['Diat,dMrtDiat','Green,dMrtGreen','DiatS1,dMrtDiatS1'],
                'DiatS1,dBurS1Diat' : ['DiatS1,dBurS1Diat'], # burial of benthic algae is currently a true sink in the model
                'ZERO: DiatS1,dResS1Diat' : ['DiatS1,dResS1Diat'], 
                'ZERO: DiatS1,dSWBuS1Dia' : ['DiatS1,dSWBuS1Dia'],
                'ZERO: DiatS1,dDigS1Diat' : ['DiatS1,dDigS1Diat'] },
     'Zoopl' : { 'Zoopl_VER_dZ_VERmor' : ['Zoopl_V,dZ_Vmor','Zoopl_E,dZ_Emor','Zoopl_R,dZ_Rmor'],
                 'Zoopl_VR_VRgr' : ['Zoopl_V,dZ_Vgr', 'Zoopl_R,dZ_Rgr']},
     'DIN' : {
             'NH4,dMinDetN' : ['NH4,dMinDetNS1', 'NH4,dMinDetNS2'], # this is a source
             'DIN,dDINUpt'  : ['NO3,dNO3Upt', 'NH4,dNH4Upt', # uptake is a sink, uptake for diatoms and greens is lumped together
                               'NH4,dNH4UptS1', 'NH4,dNH4US1D', 'NO3,dNO3UptS1'], # there are three additional uptake terms for benthic algae       
             'NO3,dDenit'   : ['NO3,dDenitWat', 'NO3,dDenitSed'], # this is a sink
             'NO3,dNiDen'     : ['NO3,dNiDen'],    # this is a sink that is usually zero but sometimes has very large spikes
             'NH4,dMinPON1' : ['NH4,dMinPON1'], # this is a source
             'NH4,dMinDON' : ['NH4,dMinDON'], # this is a source
             'NH4,dZ_NRes' : ['NH4,dZ_NRes'], # this is a source
             'NH4,dNH4Aut'  : ['NH4,dNH4Aut'], # this is a source
             'ZERO: NH4,dNH4AUTS1' : ['NH4,dNH4AUTS1'], # this is zero for now
             'ZERO: NH4,dNitrif + NO3,dNitrif'  : ['NO3,dNitrif', 'NH4,dNitrif'],     # these should cancel out    
             },
    'TN' : { 'NO3,dDenit' : ['NO3,dDenitWat', 'NO3,dDenitSed'], # this should be a SINK for TN
             'NO3,dNiDen'     : ['NO3,dNiDen'],    # this is a very small sink
             'Algae,dSedAlgae' : ['Diat,dSedDiat', 'Green,dSedGreen'],  # this should be a SINK for TN
             'PON1,dSedPON1' : ['PON1,dSedPON1'], # this should be a SINK for TN 
             'DiatS1,dMrtDiatS1' : ['DiatS1,dMrtDiatS1'], # dead benthic algae turn into detritus through the 'DetNS1,dMrtDetNS1' term, so this is a water column sink
             'DiatS1,dBurS1Diat' : ['DiatS1,dBurS1Diat'], # burial of benthic algae appears to be a true sink -- it leaves the model completely!
             'NH4,dMinDetNS' : ['NH4,dMinDetNS1', 'NH4,dMinDetNS2'], # this should be a SOURCE for TN
             'NH4,dClam_NRes' : ['NH4,dM_NRes','NH4,dG4_NRes'], # this is a SOURCE for TN (clam pee)
             'PON1,dClam_NDef' : ['PON1,dM_NDef','PON1,dG4_NDef'], # this is a SOURCE for TN (clam poo)
             'Algae,dClam_Algae' : ['Diat,dM_Diat','Diat,dG4_Diat','Green,dM_Green','Green,dG4_Green'], # this is a SINK for TN (clams eat algae)
             'PON1,dClam_PON1' : ['PON1,dM_PON1','PON1,dG4_PON1'], # this is a SINK for TN (clams eat dead algae?)
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
             'ZERO: NH4,dNH4AUTS1 + DiatS1,dResS1Diat + DiatS1,dSWBuS1Dia + DiatS1,dDigS1Diat' : ['NH4,dNH4AUTS1',      # new benthic algae terms, each one is zero, for now
                                                                                      'DiatS1,dResS1Diat', 
                                                                                      'DiatS1,dSWBuS1Dia', 
                                                                                      'DiatS1,dDigS1Diat'], 
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt' : 
                                                                                     ['Diat,dPPDiat',
                                                                                      'Diat,dcPPDiat',
                                                                                      'Green,dPPGreen',
                                                                                      'Green,dcPPGreen',
                                                                                      'NH4,dNH4Upt',
                                                                                      'NO3,dNO3Upt'],
             # identify groups of terms that sum to zero
             'ZERO: NH4,dNH4UptS1 + NH4,dNH4US1D + NO3,dNO3UptS1 + DiatS1,dPPDiatS1' : ['NH4,dNH4UptS1', # these terms sum to zero, benthic algae productivity
                                                                                        'NH4,dNH4US1D', 
                                                                                        'NO3,dNO3UptS1', 
                                                                                        'DiatS1,dPPDiatS1'], 
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
             'NO3,dNiDen'     : ['NO3,dNiDen'],    # this is a very small sink
             'DetNS2,dBurS2DetN' : ['DetNS2,dBurS2DetN'], # this appears to act as a true SINK for TN as well
             'DiatS1,dBurS1Diat' : ['DiatS1,dBurS1Diat'], # burial of benthic algae also seems to be a true sink
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
             'ZERO: NH4,dNH4AUTS1 + DiatS1,dResS1Diat + DiatS1,dSWBuS1Dia + DiatS1,dDigS1Diat' : ['NH4,dNH4AUTS1',      # new benthic algae terms, each one is zero
                                                                                                  'DiatS1,dResS1Diat', 
                                                                                                  'DiatS1,dSWBuS1Dia', 
                                                                                                  'DiatS1,dDigS1Diat'], 
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'ZERO: Diat,dPPDiat + Diat,dcPPDiat + Green,dPPGreen + Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt' : 
                                                                                     ['Diat,dPPDiat',
                                                                                      'Diat,dcPPDiat',
                                                                                      'Green,dPPGreen',
                                                                                      'Green,dcPPGreen',
                                                                                      'NH4,dNH4Upt',
                                                                                      'NO3,dNO3Upt'],
             # identify groups of terms that sum to zero
             'ZERO: NH4,dNH4UptS1 + NH4,dNH4US1D + NO3,dNO3UptS1 + DiatS1,dPPDiatS1' : ['NH4,dNH4UptS1', # these terms sum to zero, benthic algae productivity
                                                                                        'NH4,dNH4US1D', 
                                                                                        'NO3,dNO3UptS1', 
                                                                                        'DiatS1,dPPDiatS1'], 
             'ZERO: NO3,dNitrif + NH4,dNitrif' : ['NO3,dNitrif','NH4,dNitrif'],      # these two sum to zero
             'ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1' : ['DON,dCnvDPON1','PON1,dCnvDPON1'], # these two sum to zero
             'ZERO: NH4,dMinPON1 + PON1,dMinPON1' : ['NH4,dMinPON1','PON1,dMinPON1'],   # these two sum to zero
             'ZERO: NH4,dMinDON + DON,dMinDON' : ['NH4,dMinDON','DON,dMinDON'],       # these two sum to zero
             'ZERO: DetNS1,dMrtDetNS1 + DiatS1,dMrtDiatS1' : ['DetNS1,dMrtDetNS1', 'DiatS1,dMrtDiatS1'], # these new benthic algae terms sum to zero -- note currently autolysis is zero, may need to add here if it isn't later on
             'ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut' :  [ # these still sum to zero with benthic algae
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
             'ZERO: PON1,dSedPON1 + DetNS1,dSedPON1' : ['PON1,dSedPON1', # these still sum to zero with benthic algae
                                                        'DetNS1,dSedPON1'],
             'ZERO: NH4,dMinDetNS1 + NH4,dMinDetNS2 + DetNS1,dMinDetNS1 + DetNS2,dMinDetNS2' : ['NH4,dMinDetNS1',    # these still sum to zero with benthic algae
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
             'ZERO: DetNS1,dBurS1DetN + DetNS2,dBurS1DetN' : ['DetNS1,dBurS1DetN', 'DetNS2,dBurS1DetN'],  # these still sum to zero with benthic algae
             'ZERO: Diat,dSedDiat + Green,dSedGreen + DetNS1,dSedAlgN' :  ['Diat,dSedDiat','Green,dSedGreen','DetNS1,dSedAlgN'] # these still sum to zero with benthic algae
    },

    # have not gotten to these yet, leaving reactions ungrouped...
    'TotalDetNS' : {},
 #   'TP' : {},
 #   'TP_include_sediment' : {},
 #   'TotalDetPS' : {},
 #   'TotalDetSi' : {},
 #   'Grazer4' : {},
 #   'Mussel' : {},
 #   'Clams' : {}
}

# which parameters to check for mass conservation? put them in a list
# (this is used for step4_check_mass_conservation.py)
mass_cons_check_param_list = ['TN_include_sediment']

# for each parameter in mass_cons_check_param_list, pick a reaction to normalize everything by -- error in the terms
# that should be zero will be measured as a percentage of this term (this is used for step4_check_mass_conservation.py)
mass_cons_check_normalize_by_dict = {'TN_include_sediment' : 'NO3,dDenit'}

##############################
# import stuff
##############################

import sys, os
import logging
import datetime
import xarray as xr
if not stompy_dir in sys.path:
    sys.path.insert(0,stompy_dir)
import stompy.model.delft.io as dio
import socket

##############################
# functions
##############################

def get_shapefile_paths(model_inout_dir, runid, is_delta):

	'''
	poly_path, tran_path = get_shapefile_paths(model_inout_dir, runid, is_delta)

	given the path to the base directory for all our model input, the runid, and boolean saying whether or not this is a delta run,
	automatically find the path to the shapefiles used to set up the run (note we are using an earlier version of the open bay FR shapefiles, 
	which match the later version except that the later verison has some more polygons that we aren't analyzing, so saves space to leave them out
	'''

	if is_delta:

		# path to shapefiles used in final delta runs
	    tran_path =  os.path.join(model_inout_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Transects_Dave_Plus_WB_v4.shp') 
	    poly_path = os.path.join(model_inout_dir,'Delta','inputs','shapefiles','control_volumes','Delta_Fullres_Polygons_Dave_Plus_WB_v4.shp')
	
	elif 'FR' in runid:

	    ## path to the full res shapefile
	    shpfn_tran =  os.path.join(model_inout_dir,'inputs','shapefiles','Agg_exchange_lines_plus_subembayments_shoal_channel.shp') 
	    shpfn_poly = os.path.join(model_inout_dir,'inputs','shapefiles','Agg_mod_contiguous_plus_subembayments_shoal_channel.shp')

	elif 'G141' in runid:
	    
		# path to shapefiles used in aggregated grid runs
	    tran_path = os.path.join(model_inout_dir,'inputs','shapefiles','Agg_exchange_lines_141.shp')
	    poly_path =  os.path.join(model_inout_dir,'inputs','shapefiles','Agg_mod_contiguous_141.shp')
	
	return poly_path, tran_path

def get_group_def_path(runid):

	''' 
	group_def_path, group_con_path = get_group_def_path(runid)
	
	get the filenames for the group definition file and the group connectivity file, assuming they are in 
	the same repository as this script... later we will make a copy of these files into the run folder for record
	keeping purposes
	'''

	# directory where the aggregated groups are defined (this is an input directory, probably don't want to change it)
	group_def_dir = os.path.join('..','..','Definitions','group_definitions')

	# filenames
	if 'FR' in runid:
		res = 'FR'
	elif 'G141' in runid:
	    res = '141'
	else:
		raise Exception('runid %s does not fit expected pattern' % runid)

	group_def_path = os.path.abspath(os.path.join(group_def_dir,'control_volume_definitions_%s.txt' % res))
	group_con_path = os.path.abspath(os.path.join(group_def_dir,'connectivity_definitions_%s.txt' % res))

	return group_def_path, group_con_path


def get_water_year(runid):

	'''
	water_year = get_water_year(runid)

	given the runid, get the string for the water year folder the run is stored in
	'''

	if 'FR' in runid:
	    # this extracts the 2 digit water year, assuming format of runid is like FR13_003 for WY2013 run 003
	    yr = int(runid.split('_')[0][2:])
	    # turn into water year string
	    water_year = 'WY%d' % (2000 + yr)
	elif 'G141' in runid:
	    # there are two formats for agg runs, G141_13_003 is water year 2013, G141_13to18_207 is water years 2013-2018
	    # get the string that represents the water year
	    yr = runid.split('_')[1]
	    # if it spans mutlple water years (such as '13to18'), keep the string, just add 'WY' in front of it
	    if 'to' in yr:
	        water_year = 'WY' + yr
	    # otherwise extract the integer and add it to 2000
	    else:
	        water_year = 'WY%d' % (2000 + int(yr))
	else:
		raise Exception('runid %s does not fit expected pattern' % runid)

	return water_year

def get_run_dir(model_inout_dir, runid, is_delta):

	'''
	run_dir = get_run_dir(model_inout_dir, runid, is_delta)

	given the path to the base directory for all our model input, the runid, and boolean saying whether or not this is a delta run,
	automatically find the path to the directory for the run, where all the input and output files are stored
	'''

	# paths to input shapefiles, *.his and *-bal.his files, *.lsp file
	if is_delta:
	    run_dir = os.path.join(model_inout_dir,'Delta','BGC_model','Full_res',runid)
	elif 'FR' in runid:
		water_year = get_water_year(runid)
		run_dir = os.path.join(model_inout_dir,'Full_res',water_year,runid)
	elif 'G141' in runid:
		water_year = get_water_year(runid)
		run_dir = os.path.join(model_inout_dir,'Grid141',water_year,runid)
	else:
		raise Exception('runid %s does not fit expected pattern' % runid)

	return run_dir

def get_lsp_path(run_dir, is_delta):

	'''
	lsp_path = get_lsp_path(run_dir, is_delta)

	given the path to the base directory for all our model input, the runid, and boolean saying whether or not this is a delta run,
	automatically find the path to the *.lsp file
	'''

	if is_delta:
		lsp_path = os.path.join(run_dir,'%s.lsp' % water_year.lower()) 
	elif 'FR' in runid:
	    lsp_path = os.path.join(run_dir,'sfbay_dynamo000.lsp')
	elif 'G141' in runid:
	    lsp_path = os.path.join(run_dir,'sfbay_dynamo000.lsp')
	else:
		raise Exception('runid %s does not fit expected pattern' % runid)

	return lsp_path

def get_balance_table_dir(run_dir, is_delta):

	''' 
	balance_table_dir = get_balance_table_dir(run_dir, is_delta)

	unless user says otherwise, put the balance tables in the run directory, in a subfolder called "Balance_Tables"
	'''

	balance_table_dir = os.path.join(run_dir,'Balance_Tables')
	return balance_table_dir

def logger_cleanup():

    ''' check for open log files, close them, release handlers'''

    # clean up logging
    logger = logging.getLogger()
    handlers = list(logger.handlers)
    if len(handlers)>0:
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler) 

######################
# main
######################

# get the paths to the shapefiles used in the model input
if poly_path is None or tran_path is None:
	poly_path, tran_path = get_shapefile_paths(model_inout_dir, runid, is_delta)

# get the paths to the group and group connectivity definition files for the aggregated groups (used in step5_)
if group_def_path is None or group_con_path is None:
	group_def_path, group_con_path = get_group_def_path(runid)

# get the run directory
if run_dir is None:
	run_dir = get_run_dir(model_inout_dir, runid, is_delta)

# get the path to the lsp file
if lsp_path is None:
	lsp_path = get_lsp_path(run_dir, is_delta)

# get the balance table directory
if balance_table_dir is None:
	balance_table_dir = get_balance_table_dir(run_dir, is_delta)

# create the balance table directory if it does not exist
if not os.path.exists(balance_table_dir):
	print('Creating folder to store balance tables here:')
	print('    %s' % balance_table_dir)
	os.makedirs(balance_table_dir)
else: 
	print('Balance table folder already exists here:')
	print('    %s' % balance_table_dir)

# setup logging to file and to screen 
logger_cleanup()
logging.basicConfig(
level=logging.INFO,
format="%(asctime)s [%(levelname)s] %(message)s",
handlers=[
    logging.FileHandler(os.path.join(balance_table_dir,"log_step0.log"),'w'),
    logging.StreamHandler(sys.stdout)
])

# add some basic info to log file
user = os.getlogin()
hostname = socket.gethostname()
scriptname= __file__
conda_env=os.environ['CONDA_DEFAULT_ENV']
today= datetime.datetime.now().strftime('%b %d, %Y')
logging.info('Balance table configuration was specified on %s by %s on %s in %s using %s' % (today, user, hostname, conda_env, scriptname))

	
# print the values of all the parameters we are going to use to create the balance tables
logging.info('')
logging.info('Setting things up in step0_config.py ...')
logging.info('')
logging.info('  Processing the following run:')
logging.info('      step0_config.run_dir = %s' % run_dir)
logging.info('  Balance tables will be saved here:')
logging.info('      step0_config.balance_table_dir = %s' % balance_table_dir)
logging.info('  Delete any balance tables found in this directory before proceeding?')
logging.info('      step0_config.delete_balance_tables = %s' % delete_balance_tables)
logging.info('  Processing the following modeled substances:')
logging.info('      step0_config.substance_list = %s' % substance_list)
logging.info('  Assuming polygons and transects in the *.his and *-bal.his files are based on these shapefiles:')
logging.info('      step0_config.poly_path = %s' % poly_path)
logging.info('      step0_config.tran_path = %s' % tran_path)
logging.info('  Using the following group and connectivity definitions for aggregating in space:')
logging.info('      step0_config.group_def_path = %s' % group_def_path)
logging.info('      step0_config.group_con_path = %s' % group_con_path)
logging.info('  Using the *.lsp file here:')
logging.info('      step0_config.lsp_path = %s' % lsp_path)
logging.info('  Using the following precision in balance tables:')
logging.info('      step0_config.float_format = %s' % float_format)
logging.info('  Abort for mass conservation error?')
logging.info('      step0_config.abort_for_mass_cons_error = %s' % abort_for_mass_cons_error)
logging.info('  Tolerance for mass conservation error:')
logging.info('      step0_config.error_tol_percent = %f percent' % error_tol_percent)
logging.info('')
logger_cleanup()

