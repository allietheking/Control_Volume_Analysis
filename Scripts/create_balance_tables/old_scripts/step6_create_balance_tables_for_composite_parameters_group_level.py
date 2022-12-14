
'''
Starting with raw parameters found in 'Balance_Table_All_Parameters_By_Group.csv' file, 
calculates fluxes and rates for composite parameters defined in the user input section by
constituents and reaction groupings.
'''


#################################################################
# USER INPUT
#################################################################

# directory containting balance tables, including path. note if there are multiple runs,
# there is another layer of run folders inside this directory before the actual balance tables
#run_folders = ['G141_13to18_125'] 
#run_folders = ['FR13_021'] 
#run_folders = ['FR17_014'] 
#run_folders = ['FR18_004'] 
#run_folders = ['FR13_003'] 
run_folders = ['FR17_003']

# output directory, including path
output_dir = '/richmondvol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/'

# does run have greens?
hasgreens = False

if hasgreens:

    # composite parameters -- for each parameter, include a stoiciometric multiplier;
    # for example for diatoms, the fraction contributing to nitrogen is 0.14,
    # so when computing total nitrogen, need to multiply diatom masses by 0.14
    composite_parameters = {
     'TN'  : [[1,'NH4'], [1,'NO3'], [1,'DON'], [1,'PON1'], [0.16,'Diat'], [0.16,'Green'], [0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
     'DIN' : [[1.00, 'NH4'], [1.00, 'NO3']],
     'PON' : [[1.00, 'PON1']],
     'N-Algae' : [[0.15, 'Diat'], [0.16, 'Green']],
     'N-Zoopl' : [[0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
     'DetN' : [[1.00,'DetNS1'], [1.00,'DetNS2']],
     'Zoopl' : [[1.0,'Zoopl_V'], [1.0,'Zoopl_R'], [1.0,'Zoopl_E']],
     'TP' : [[1,'PO4'], [1,'DOP'], [1,'POP1'], [0.01,'Diat'], [0.01, 'Green'], [0.0263,'Zoopl_V'], [0.0263,'Zoopl_R'], [0.0263,'Zoopl_E']],
     'DO' : [[1,'OXY']]
    }

else:

    # composite parameters -- for each parameter, include a stoiciometric multiplier;
    # for example for diatoms, the fraction contributing to nitrogen is 0.14,
    # so when computing total nitrogen, need to multiply diatom masses by 0.14
    composite_parameters = {
     'TN'  : [[1,'NH4'], [1,'NO3'], [1,'DON'], [1,'PON1'], [0.15,'Diat'], [0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
     'DIN' : [[1.00, 'NH4'], [1.00, 'NO3']],
     'PON' : [[1.00, 'PON1']],
     'N-Algae' : [[0.15, 'Diat']],
     'N-Zoopl' : [[0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
     'DetN' : [[1.00,'DetNS1'], [1.00,'DetNS2']],
     'Zoopl' : [[1.0,'Zoopl_V'], [1.0,'Zoopl_R'], [1.0,'Zoopl_E']],
     'TP' : [[1,'PO4'], [1,'DOP'], [1,'POP1'], [0.01,'Diat'], [0.0263,'Zoopl_V'], [0.0263,'Zoopl_R'], [0.0263,'Zoopl_E']],
     'DO' : [[1,'OXY']]
    }


if hasgreens:
          
    # composite reactions -- this script will add the list of terms on the right of the colon
    # to create the term on the left of the colon. all other terms will be left alone. note the
    # plotting script drops terms that are zero, so adding terms that add to zero will create
    # a composite term that is equal to zero, and thus does not get plotted, provided the terms
    # really do sum to zero
    composite_reactions = {
   'TN' : { # the sources and sinks will be included automatically, and if you explicitly call them out WITHOUT changing the 
             # name the code will screw up and delete that term, so leave out the sources and sinks, but list them as comments here,
             # and it is ok to group some of them provided we give them NEW names that don't match the names of the original terms
             'NO3,dDenit (Mg/d)' : ['NO3,dDenitWat (Mg/d)','NO3,dDenitSed (Mg/d)'], # this should be a SINK for TN
             'N-Algae,dSedDiat (Mg/d)' : ['Diat,dSedDiat (Mg/d)', 'Green,dSedGreen (Mg/d)'], # this should be a SINK for TN
             #'PON1,dSedPON1 (Mg/d)' : ['PON1,dSedPON1 (Mg/d)'], # this should be a SINK for TN -- do not include it b/c don't want to rename it
             'NH4,dMinDetNS (Mg/d)' : ['NH4,dMinDetNS1 (Mg/d)', 'NH4,dMinDetNS2 (Mg/d)'], # this should be a SOURCE for TN
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'SHOULD SUM TO ZERO: N-Diat,dPPDiat + N-Diat,dcPPDiat + N-Green,dPPGreen + N-Green,dcPPGreen + NH4,dNH4Upt + NO3,dNO3Upt (Mg/d)' : ['Diat,dPPDiat (Mg/d)',
                                                                                      'Diat,dcPPDiat (Mg/d)',
                                                                                      'Green,dPPGreen (Mg/d)',
                                                                                      'Green,dcPPGreen (Mg/d)',
                                                                                      'NH4,dNH4Upt (Mg/d)',
                                                                                      'NO3,dNO3Upt (Mg/d)'],
             # identify groups of terms that sum to zero
             'SHOULD SUM TO ZERO: NO3,dNitrif + NH4,dNitrif (Mg/d)' : ['NO3,dNitrif (Mg/d)','NH4,dNitrif (Mg/d)'],      # these two sum to zero
             'SHOULD SUM TO ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1 (Mg/d)' : ['DON,dCnvDPON1 (Mg/d)','PON1,dCnvDPON1 (Mg/d)'], # these two sum to zero
             'SHOULD SUM TO ZERO: NH4,dMinPON1 + PON1,dMinPON1 (Mg/d)' : ['NH4,dMinPON1 (Mg/d)','PON1,dMinPON1 (Mg/d)'],   # these two sum to zero
             'SHOULD SUM TO ZERO: NH4,dMinDON + DON,dMinDON (Mg/d)' : ['NH4,dMinDON (Mg/d)','DON,dMinDON (Mg/d)'],       # these two sum to zero
             'THIS IS ZERO: NO3,dNiDen (Mg/d)'     : ['NO3,dNiDen (Mg/d)'],                             # this is zero for the current model (because no sediment model??? not sure why)
             'EACH ONE IS ZERO: PON1,dCnvPPON1 + PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN + PON1,dZ_NMrt + PON1,dZ_PON1 (Mg/d)'   : [
                                          'PON1,dCnvPPON1 (Mg/d)',  
                                          'PON1,dResS1DetN (Mg/d)', 
                                          'PON1,dResS2DetN (Mg/d)', 
                                          'PON1,dResS1DiDN (Mg/d)', 
                                          'PON1,dResS2DiDN (Mg/d)', 
                                          'PON1,dZ_NMrt (Mg/d)',      
                                          'PON1,dZ_PON1 (Mg/d)'],       
             'EACH ONE IS ZERO: Zoopl_V,dZ_Vmor, Zoopl_R,dZ_Rmor (Mg/d)' : [                                        
                                         'Zoopl_V,dZ_Vmor (Mg/d)',  
                                          'Zoopl_R,dZ_Rmor (Mg/d)',  
                                          'Zoopl_E,dZ_Emor (Mg/d)'], 
             'SHOULD SUM TO ZERO: Diat,dMrtDiat + Green,dMrtGreen + PON1,dMortDetN + NH4,dNH4Aut (Mg/d)' :  [ 
                                                'Diat,dMrtDiat (Mg/d)', 
                                                'Green,dMrtGreen (Mg/d)',
                                                'PON1,dMortDetN (Mg/d)',
                                                'NH4,dNH4Aut (Mg/d)'], 
             'SHOULD SUM TO ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Green,dZ_Grn + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec (Mg/d)' : [
                                         'NH4,dZ_NRes (Mg/d)',
                                         'PON1,dZ_NDef (Mg/d)',
                                         'PON1,dZ_NSpDet (Mg/d)',
                                         'Diat,dZ_Diat (Mg/d)',
                                         'Green,dZ_Grn (Mg/d)',
                                         'Zoopl_V,dZ_Vgr (Mg/d)',
                                         'Zoopl_R,dZ_SpwDet (Mg/d)',
                                         'Zoopl_R,dZ_Rgr (Mg/d)',
                                         'Zoopl_E,dZ_Ea (Mg/d)',
                                         'Zoopl_E,dZ_Ec (Mg/d)'
                                         ]
          },
    'DIN' : {
             'NH4,dMinDetN (Mg/d)' : ['NH4,dMinDetNS1 (Mg/d)', 'NH4,dMinDetNS2 (Mg/d)'],
             'DIN,dDINUpt (Mg/d)'  : ['NO3,dNO3Upt (Mg/d)', 'NH4,dNH4Upt (Mg/d)'],       
             'DIN,dNitrif (Mg/d)'  : ['NO3,dNitrif (Mg/d)', 'NH4,dNitrif (Mg/d)']       # these should cancel out
            },
    'PON' : {},
    'N-Algae' : {
                'N-Algae,dPPAlgae + N-Algae,dcPPAlgae (Mg/d)' : ['Diat,dPPDiat (Mg/d)','Diat,dcPPDiat (Mg/d)','Green,dPPGreen (Mg/d)','Green,dcPPGreen (Mg/d)'],
                'N-Algae,dMrtAlgae (Mg/d)' : ['Diat,dMrtDiat (Mg/d)', 'Green,dMrtGreen (Mg/d)'],
                'N-Algae,dSedAlgae (Mg/d)' : ['Diat,dSedDiat (Mg/d)', 'Green,dSedGreen (Mg/d)'],
                'N-Algae,dZ_Algae (Mg/d)'  : ['Diat,dZ_Diat (Mg/d)', 'Green,dZ_Grn (Mg/d)']
               },
    'N-Zoopl' : {'N-Zoopl_E,dZ_Ea (Mg/d)' : ['Zoopl_E,dZ_Ea (Mg/d)'],
                 'N-Zoopl_E,dZ_Ec (Mg/d)' : ['Zoopl_E,dZ_Ec (Mg/d)'],
                 'N-Zoopl_R,dZ_SpwDet (Mg/d)' : ['Zoopl_R,dZ_SpwDet (Mg/d)'],
                 'N-Zoopl_R,dZ_Rgr + N-Zoopl_V,dZ_Vgr (Mg/d)' : ['Zoopl_R,dZ_Rgr (Mg/d)','Zoopl_V,dZ_Vgr (Mg/d)']
                 },
    'DetN' :    {'DetN,dBurDetN (Mg/d)' : ['DetNS1,dBurS1DetN (Mg/d)', 'DetNS2,dBurS1DetN (Mg/d)'],
                 'DetN,dDigDetN (Mg/d)' : ['DetNS1,dDigS1DetN (Mg/d)', 'DetNS2,dDigS1DetN (Mg/d)']},
    'Zoopl' : {'Zoopl_R,dZ_Rgr + Zoopl_V,dZ_Vgr (Mg/d)' : ['Zoopl_R,dZ_Rgr (Mg/d)','Zoopl_V,dZ_Vgr (Mg/d)']},
    'TP'    : {'P-Algae,dSedAlgae (Mg/d)' : ['Diat,dSedDiat (Mg/d)','Green,dSedGreen (Mg/d)'],
               'TP-Internal (Mg/d)'     : ['Diat,dZ_Diat (Mg/d)','Green,dZ_Grn (Mg/d)',
                                         #'Diat,dSedDiat (Mg/d)', # flagged by exploratory algorithm in v1 for TN
                                         'Diat,dPPDiat (Mg/d)', 'Green,dPPGreen (Mg/d)',# part of uptake balance
                                         'Diat,dMrtDiat (Mg/d)','Green,dMrtGreen (Mg/d)',
                                         'Diat,dcPPDiat (Mg/d)','Green,dcPPGreen (Mg/d)', # part of uptake balance
                                         'Zoopl_V,dZ_Vgr (Mg/d)',
                                         'Zoopl_V,dZ_Vmor (Mg/d)',
                                         'Zoopl_R,dZ_SpwDet (Mg/d)',
                                         'Zoopl_R,dZ_Rgr (Mg/d)',
                                         'Zoopl_R,dZ_Rmor (Mg/d)',
                                         'Zoopl_E,dZ_Ea (Mg/d)',
                                         'Zoopl_E,dZ_Ec (Mg/d)',
                                         'Zoopl_E,dZ_Emor (Mg/d)',
                                         "PO4,dMinPOP1 (Mg/d)",
                                         "PO4,dMinDOP (Mg/d)",
                                         "PO4,dZ_PRes (Mg/d)",
                                         #"PO4,dMinDetPS1 (Mg/d)", # should be a source for TP
                                         #"PO4,dMinDetPS2 (Mg/d)", # should be a source for TP
                                         "PO4,dPO4Aut (Mg/d)",
                                         "PO4,dPO4Upt (Mg/d)", # part of uptake balance
                                         "POP1,dCnvPPOP1 (Mg/d)",
                                         "POP1,dCnvDPOP1 (Mg/d)",
                                         "POP1,dMinPOP1 (Mg/d)",
                                         "POP1,dZ_PMrt (Mg/d)",
                                         "POP1,dZ_PDef (Mg/d)",
                                         "POP1,dZ_PSpDet (Mg/d)",
                                         "POP1,dZ_POP (Mg/d)",
                                         #"POP1,dSedPOP1 (Mg/d)", # should be a sink for TP
                                         "POP1,dMortDetP (Mg/d)",
                                         #"POP1,dResS1DetP (Mg/d)", # should be a source for TP
                                         #"POP1,dResS2DetP (Mg/d)", # should be a source for TP
                                         "POP1,dResS1DiDP (Mg/d)",
                                         "POP1,dResS2DiDP (Mg/d)",
                                         "DOP,dCnvDPOP1 (Mg/d)",
                                         "DOP,dMinDOP (Mg/d)"]
                 },
                                         
    'DO' : {}
    }

else:
          
    composite_reactions = {
    'TN' : { # the sources and sinks will be included automatically, and if you explicitly call them out WITHOUT changing the 
             # name the code will screw up and delete that term, so leave out the sources and sinks, but list them as comments here,
             # and it is ok to group some of them provided we give them NEW names that don't match the names of the original terms
             'NO3,dDenit (Mg/d)' : ['NO3,dDenitWat (Mg/d)','NO3,dDenitSed (Mg/d)'], # this should be a SINK for TN
             'N-Diat,dSedDiat (Mg/d)' : ['Diat,dSedDiat (Mg/d)'], # this should be a SINK for TN
             #'PON1,dSedPON1 (Mg/d)' : ['PON1,dSedPON1 (Mg/d)'], # this should be a SINK for TN -- do not include it b/c don't want to rename it
             'NH4,dMinDetNS (Mg/d)' : ['NH4,dMinDetNS1 (Mg/d)', 'NH4,dMinDetNS2 (Mg/d)'], # this should be a SOURCE for TN
             # in the PRE August 2020 model this does NOT sum to zero but it should
             'SHOULD SUM TO ZERO: N-Diat,dPPDiat + N-Diat,dcPPDiat + NH4,dNH4Upt + NO3,dNO3Upt (Mg/d)' : ['Diat,dPPDiat (Mg/d)',
                                                                                      'Diat,dcPPDiat (Mg/d)',
                                                                                      'NH4,dNH4Upt (Mg/d)',
                                                                                      'NO3,dNO3Upt (Mg/d)'],
             # identify groups of terms that sum to zero
             'SHOULD SUM TO ZERO: NO3,dNitrif + NH4,dNitrif (Mg/d)' : ['NO3,dNitrif (Mg/d)','NH4,dNitrif (Mg/d)'],      # these two sum to zero
             'SHOULD SUM TO ZERO: DON,dCnvDPON1 + PON1,dCnvDPON1 (Mg/d)' : ['DON,dCnvDPON1 (Mg/d)','PON1,dCnvDPON1 (Mg/d)'], # these two sum to zero
             'SHOULD SUM TO ZERO: NH4,dMinPON1 + PON1,dMinPON1 (Mg/d)' : ['NH4,dMinPON1 (Mg/d)','PON1,dMinPON1 (Mg/d)'],   # these two sum to zero
             'SHOULD SUM TO ZERO: NH4,dMinDON + DON,dMinDON (Mg/d)' : ['NH4,dMinDON (Mg/d)','DON,dMinDON (Mg/d)'],       # these two sum to zero
             'THIS IS ZERO: NO3,dNiDen (Mg/d)'     : ['NO3,dNiDen (Mg/d)'],                             # this is zero for the current model (because no sediment model??? not sure why)
             'EACH ONE IS ZERO: PON1,dCnvPPON1 + PON1,dResS1DetN + PON1,dResS2DetN + PON1,dResS1DiDN + PON1,dResS2DiDN + PON1,dZ_NMrt + PON1,dZ_PON1 (Mg/d)'   : [
                                          'PON1,dCnvPPON1 (Mg/d)',  
                                          'PON1,dResS1DetN (Mg/d)', 
                                          'PON1,dResS2DetN (Mg/d)', 
                                          'PON1,dResS1DiDN (Mg/d)', 
                                          'PON1,dResS2DiDN (Mg/d)', 
                                          'PON1,dZ_NMrt (Mg/d)',      
                                          'PON1,dZ_PON1 (Mg/d)'],       
             'EACH ONE IS ZERO: Zoopl_V,dZ_Vmor, Zoopl_R,dZ_Rmor (Mg/d)' : [                                        
                                         'Zoopl_V,dZ_Vmor (Mg/d)',  
                                          'Zoopl_R,dZ_Rmor (Mg/d)',  
                                          'Zoopl_E,dZ_Emor (Mg/d)'], 
             'SHOULD SUM TO ZERO: Diat,dMrtDiat + PON1,dMortDetN + NH4,dNH4Aut (Mg/d)' :  [ 
                                                'Diat,dMrtDiat (Mg/d)',
                                                'PON1,dMortDetN (Mg/d)',
                                                'NH4,dNH4Aut (Mg/d)'],
             'SHOULD SUM TO ZERO: NH4,dZ_NRes + PON1,dZ_NDef + PON1,dZ_NSpDet + Diat,dZ_Diat + Zoopl_V,dZ_Vgr + Zoopl_R,dZ_SpwDet + Zoopl_R,dZ_Rgr + Zoopl_E,dZ_Ea + Zoopl_E,dZ_Ec (Mg/d) (Mg/d)' : [
                                         'NH4,dZ_NRes (Mg/d)',
                                         'PON1,dZ_NDef (Mg/d)',
                                         'PON1,dZ_NSpDet (Mg/d)',
                                         'Diat,dZ_Diat (Mg/d)',
                                         'Zoopl_V,dZ_Vgr (Mg/d)',
                                         'Zoopl_R,dZ_SpwDet (Mg/d)',
                                         'Zoopl_R,dZ_Rgr (Mg/d)',
                                         'Zoopl_E,dZ_Ea (Mg/d)',
                                         'Zoopl_E,dZ_Ec (Mg/d)'
                                         ]
          },
    'DIN' : {
             'NH4,dMinDetN (Mg/d)' : ['NH4,dMinDetNS1 (Mg/d)', 'NH4,dMinDetNS2 (Mg/d)'],
             'DIN,dDINUpt (Mg/d)'  : ['NO3,dNO3Upt (Mg/d)', 'NH4,dNH4Upt (Mg/d)'],       
             'DIN,dNitrif (Mg/d)'  : ['NO3,dNitrif (Mg/d)', 'NH4,dNitrif (Mg/d)']       # these should cancel out
            },
    'PON' : {},
    'N-Algae' : {
                'N-Diat,dPPDiat + N-Diat,dcPPDiat (Mg/d)' : ['Diat,dPPDiat (Mg/d)','Diat,dcPPDiat (Mg/d)'],
                'N-Diat,dMrtDiat (Mg/d)' : ['Diat,dMrtDiat (Mg/d)'],
                'N-Diat,dSedDiat (Mg/d)' : ['Diat,dSedDiat (Mg/d)'],
                'N-Diat,dZ_Diat (Mg/d)'  : ['Diat,dZ_Diat (Mg/d)']
               },
    'N-Zoopl' : {'N-Zoopl_E,dZ_Ea (Mg/d)' : ['Zoopl_E,dZ_Ea (Mg/d)'],
                 'N-Zoopl_E,dZ_Ec (Mg/d)' : ['Zoopl_E,dZ_Ec (Mg/d)'],
                 'N-Zoopl_R,dZ_SpwDet (Mg/d)' : ['Zoopl_R,dZ_SpwDet (Mg/d)'],
                 'N-Zoopl_R,dZ_Rgr + N-Zoopl_V,dZ_Vgr (Mg/d)' : ['Zoopl_R,dZ_Rgr (Mg/d)','Zoopl_V,dZ_Vgr (Mg/d)']
                 },
    'DetN' : {},
    'Zoopl' : {'Zoopl_R,dZ_Rgr + Zoopl_V,dZ_Vgr (Mg/d)' : ['Zoopl_R,dZ_Rgr (Mg/d)','Zoopl_V,dZ_Vgr (Mg/d)']},
    'TP'    : {'P-Diat,dSedDiat (Mg/d)' : ['Diat,dSedDiat (Mg/d)'],
               'TP-Internal (Mg/d)'     : ['Diat,dZ_Diat (Mg/d)',
                                         #'Diat,dSedDiat (Mg/d)', # flagged by exploratory algorithm in v1 for TN
                                         'Diat,dPPDiat (Mg/d)', # part of uptake balance
                                         'Diat,dMrtDiat (Mg/d)',
                                         'Diat,dcPPDiat (Mg/d)', # part of uptake balance
                                         'Zoopl_V,dZ_Vgr (Mg/d)',
                                         'Zoopl_V,dZ_Vmor (Mg/d)',
                                         'Zoopl_R,dZ_SpwDet (Mg/d)',
                                         'Zoopl_R,dZ_Rgr (Mg/d)',
                                         'Zoopl_R,dZ_Rmor (Mg/d)',
                                         'Zoopl_E,dZ_Ea (Mg/d)',
                                         'Zoopl_E,dZ_Ec (Mg/d)',
                                         'Zoopl_E,dZ_Emor (Mg/d)',
                                         "PO4,dMinPOP1 (Mg/d)",
                                         "PO4,dMinDOP (Mg/d)",
                                         "PO4,dZ_PRes (Mg/d)",
                                         #"PO4,dMinDetPS1 (Mg/d)", # should be a source for TP
                                         #"PO4,dMinDetPS2 (Mg/d)", # should be a source for TP
                                         "PO4,dPO4Aut (Mg/d)",
                                         "PO4,dPO4Upt (Mg/d)", # part of uptake balance
                                         "POP1,dCnvPPOP1 (Mg/d)",
                                         "POP1,dCnvDPOP1 (Mg/d)",
                                         "POP1,dMinPOP1 (Mg/d)",
                                         "POP1,dZ_PMrt (Mg/d)",
                                         "POP1,dZ_PDef (Mg/d)",
                                         "POP1,dZ_PSpDet (Mg/d)",
                                         "POP1,dZ_POP (Mg/d)",
                                         #"POP1,dSedPOP1 (Mg/d)", # should be a sink for TP
                                         "POP1,dMortDetP (Mg/d)",
                                         #"POP1,dResS1DetP (Mg/d)", # should be a source for TP
                                         #"POP1,dResS2DetP (Mg/d)", # should be a source for TP
                                         "POP1,dResS1DiDP (Mg/d)",
                                         "POP1,dResS2DiDP (Mg/d)",
                                         "DOP,dCnvDPOP1 (Mg/d)",
                                         "DOP,dMinDOP (Mg/d)"]
                 },
                                         
    'DO' : {}
    }
  
#################################################################
# IMPORT PACKAGES
#################################################################

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import datetime as dt

#########################
# MAIN
#########################

# loop through run folders
for run_folder in run_folders:

    # input file -- balance table for raw parameters, aggregated into groups
    input_balance_table_filename = os.path.join(output_dir,run_folder,'Balance_Table_All_Parameters_By_Group.csv')
    
    # prefix for output balance tables -- balance tables for composite parameters will be saved individually
    output_balance_table_prefix = os.path.join(output_dir,run_folder,'Balance_Table_By_Group_Composite_Parameter_')
    
    # read main balance table with raw parameters
    df_master = pd.read_csv(input_balance_table_filename)
    print(set(df_master.group))
    # loop through composite parameters
    for cp in composite_parameters.keys():
    
        # get list of constituent parameters and corresponding multipliers
        mult_list = []
        param_list = []
        for mp in composite_parameters[cp]:
            mult_list.append(mp[0])
            param_list.append(mp[1])
    
        # make a list of columns that includes the groups, time stamps, volume
        # and area information, and all the terms corresponding to the constituent
        # parameters for this composite parameter
        columns = ['group', 'time', 'Volume (m^3)', 'Area (m^2)']
        for col in df_master.columns:
            if ',' in col:
                param = col.split(',')[0]
                if param in param_list:
                    columns.append(col)
                    
        # make a data frame from the subset of columns identified above
        df = df_master[columns].copy(deep=True)
        
        # loop through the columns and apply the multiplier for each parameter
        for mult, param in zip(mult_list, param_list):
            for col in df.columns:
                if ',' in col:
                    if param == col.split(',')[0]:
                        df[col] = df[col].values.copy() * mult
        
        # insert new columns for the composite parameter mass budget terms, all 
        # except the reactions, initializing values at zero
        df.insert(4,cp + ',Mass (Mg)',0.0)
        df.insert(5,cp + ',dMass/dt (Mg/d)',0.0)
        df.insert(6,cp + ',Net Flux In (Mg/d)',0.0)
        df.insert(7,cp + ',Net Load (Mg/d)',0.0)
        df.insert(8,cp + ',Net Transport In (Mg/d)',0.0)
        df.insert(9,cp + ',Net Reaction (Mg/d)',0.0)
        df.insert(10,cp + ',dMass/dt, Balance Check (Mg/d)',0.0)
        df.insert(11,cp + ',Flux In from N (Mg/d)',0.0)
        df.insert(12,cp + ',Flux In from S (Mg/d)',0.0) 
        df.insert(13,cp + ',Flux In from E (Mg/d)',0.0)
        df.insert(14,cp + ',Flux In from W (Mg/d)',0.0)
        
        # loop through the composite parameter and add the contributions to the mass budget terms
        for param in param_list:
        
            # add contribution from each constituent parameter to the composite parameter
            df[cp + ',Mass (Mg)'] += df[param + ',Mass (Mg)'].values.copy()
            df[cp + ',dMass/dt (Mg/d)'] += df[param + ',dMass/dt (Mg/d)'].values.copy()
            df[cp + ',Net Flux In (Mg/d)'] += df[param + ',Net Flux In (Mg/d)'].values.copy()
            df[cp + ',Net Load (Mg/d)'] += df[param + ',Net Load (Mg/d)'].values.copy()
            df[cp + ',Net Transport In (Mg/d)'] += df[param + ',Net Transport In (Mg/d)'].values.copy()
            df[cp + ',Net Reaction (Mg/d)'] += df[param + ',Net Reaction (Mg/d)'].values.copy()
            df[cp + ',dMass/dt, Balance Check (Mg/d)'] += df[param + ',dMass/dt, Balance Check (Mg/d)'].values.copy() 
            df[cp + ',Flux In from N (Mg/d)'] += df[param + ',Flux In from N (Mg/d)'].values.copy()
            df[cp + ',Flux In from S (Mg/d)'] += df[param + ',Flux In from S (Mg/d)'].values.copy()
            df[cp + ',Flux In from E (Mg/d)'] += df[param + ',Flux In from E (Mg/d)'].values.copy()
            df[cp + ',Flux In from W (Mg/d)'] += df[param + ',Flux In from W (Mg/d)'].values.copy()
            
            # then drop all the columns from that consituent parameter
            df.drop(columns = [param + ',Mass (Mg)', 
                               param + ',dMass/dt (Mg/d)',
                               param + ',Net Flux In (Mg/d)',
                               param + ',Net Load (Mg/d)',
                               param + ',Net Transport In (Mg/d)',
                               param + ',Net Reaction (Mg/d)',
                               param + ',dMass/dt, Balance Check (Mg/d)',
                               param + ',Flux In from N (Mg/d)',
                               param + ',Flux In from S (Mg/d)',
                               param + ',Flux In from E (Mg/d)',
                               param + ',Flux In from W (Mg/d)'], inplace=True)
       

        # loop through the composite reaction terms
        for cr in composite_reactions[cp].keys():

            # make sure none of the component reactions have the same name as the composite reaction
            if cr in composite_reactions[cp][cr]:
                print('Error! The composite reaction %s in the %s budget has the SAME NAME as' % (cr,cp))
                print('one of the component reactions that contribute to it. Here is the list of those component reactions:')
                print(composite_reactions[cp][cr])
                print('This will cause BIG PROBLEMS so please give it a different name')
                raise Exception('REACTION NAME ERROR')
        
            # add a column for this composite reaction term, initializing at zero
            df[cr] = 0.0
            
            # then loop through all the constituent reactions and add them, deleting afterwards
            for cri in composite_reactions[cp][cr]:
                df[cr] += df[cri].values.copy()
                df = df.drop(columns = [cri])
        
        # finally, save the composite parameter budget to its own data frame
        df.to_csv(output_balance_table_prefix + cp + '.csv', index=False)
        
            



