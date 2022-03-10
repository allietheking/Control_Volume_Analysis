# composite parameters -- for each parameter, include a stoiciometric multiplier;
# for example for diatoms, the fraction contributing to nitrogen is 0.14,
# so when computing total nitrogen, need to multiply diatom masses by 0.14
composite_parameters = {
 'TN'  : [[1,'NH4'], [1,'NO3'], [1,'DON'], [1,'PON1'], [0.15,'Diat'], [0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
 'DIN' : [[1.00, 'NH4'], [1.00, 'NO3']],
 'PON' : [[1.00, 'PON1']],
 'N-Diat' : [[0.15, 'Diat']],
 'N-Zoopl' : [[0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']],
 'DetN' : [[1.00,'DetNS1'], [1.00,'DetNS2']],
 'Zoopl' : [[1.0,'Zoopl_V'], [1.0,'Zoopl_R'], [1.0,'Zoopl_E']],
 'TP' : [[1,'PO4'], [1,'DOP'], [1,'POP1'], [0.01,'Diat'], [0.0263,'Zoopl_V'], [0.0263,'Zoopl_R'], [0.0263,'Zoopl_E']],
 'DO' : [[1,'OXY']]
}

          
# when developing the balance table for a composite parameter, the aggregation script by default includes every 
# reaciton term for every component parameter. we may want to group some of these component reaction terms together 
# and give the grouped term a new name, to declutter the list of reactions for the component parameter. some of these 
# sets of component reaction terms should sum to zero. here we define groups of parameters, and later the code will
# drop the terms that sum to zero ... all the component reactions not mentioned in this list will be preserved as is
composite_reactions = {
'TN' : { 'N-Diat,dPPDiat + N-Diat,dcPPDiat + NH4,dNH4Upt + NO3,dNO3Upt (Mg/d)' : ['Diat,dPPDiat (Mg/d)',    
                                                                                  'Diat,dcPPDiat (Mg/d)',
                                                                                  'NH4,dNH4Upt (Mg/d)',
                                                                                  'NO3,dNO3Upt (Mg/d)'],
         'N-Diat,dSedDiat (Mg/d)' : ['Diat,dSedDiat (Mg/d)'],
         'TN-Internal (Mg/d)'     : ['Diat,dZ_Diat (Mg/d)',                                                  
                                     #'Diat,dSedDiat (Mg/d)', # flagged by exploratory algorithm in v1
                                     #'Diat,dPPDiat (Mg/d)', # part of uptake balance
                                     'Diat,dMrtDiat (Mg/d)',
                                     #'Diat,dcPPDiat (Mg/d)', # part of uptake balance
                                     'Zoopl_V,dZ_Vgr (Mg/d)',
                                     'Zoopl_V,dZ_Vmor (Mg/d)',
                                     'Zoopl_R,dZ_SpwDet (Mg/d)',
                                     'Zoopl_R,dZ_Rgr (Mg/d)',
                                     'Zoopl_R,dZ_Rmor (Mg/d)',
                                     'Zoopl_E,dZ_Ea (Mg/d)',
                                     'Zoopl_E,dZ_Ec (Mg/d)',
                                     'Zoopl_E,dZ_Emor (Mg/d)',
                                     #'NO3,dDenitWat (Mg/d)', # this should be a SINK for TN
                                     'NO3,dNitrif (Mg/d)',
                                     #'NO3,dDenitSed (Mg/d)', # this should be a SINK for TN
                                     'NO3,dNiDen (Mg/d)',
                                     #'NO3,dNO3Upt (Mg/d)',  # part of uptake balance
                                     'NH4,dMinPON1 (Mg/d)',
                                     'NH4,dMinDON (Mg/d)',
                                     'NH4,dNitrif (Mg/d)',
                                     #'NH4,dMinDetNS1 (Mg/d)', # this should be a SOURCE for TN
                                     #'NH4,dMinDetNS2 (Mg/d)', # this should be a SOURCE for TN
                                     'NH4,dZ_NRes (Mg/d)',
                                     'NH4,dNH4Aut (Mg/d)',
                                     #'NH4,dNH4Upt (Mg/d)',    # part of uptake balance
                                     'DON,dCnvDPON1 (Mg/d)',
                                     'DON,dMinDON (Mg/d)',
                                     'PON1,dCnvPPON1 (Mg/d)',
                                     'PON1,dCnvDPON1 (Mg/d)',
                                     'PON1,dMinPON1 (Mg/d)',
                                     'PON1,dZ_NMrt (Mg/d)',
                                     'PON1,dZ_NDef (Mg/d)',
                                     'PON1,dZ_NSpDet (Mg/d)',
                                     'PON1,dZ_PON1 (Mg/d)',
                                     #'PON1,dSedPON1 (Mg/d)', # this should be a SINK for TN
                                     'PON1,dMortDetN (Mg/d)',
                                     #'PON1,dResS1DetN (Mg/d)', # this should be a SOURCE for TN
                                     #'PON1,dResS2DetN (Mg/d)', # this should be a SOURCE for TN
                                     'PON1,dResS1DiDN (Mg/d)',
                                     'PON1,dResS2DiDN (Mg/d)']
      },
'DIN' : {
         'NH4,dMinDetN (Mg/d)' : ['NH4,dMinDetNS1 (Mg/d)', 'NH4,dMinDetNS2 (Mg/d)'],
         'DIN,dDINUpt (Mg/d)'  : ['NO3,dNO3Upt (Mg/d)', 'NH4,dNH4Upt (Mg/d)'],       
         'DIN,dNitrif (Mg/d)'  : ['NO3,dNitrif (Mg/d)', 'NH4,dNitrif (Mg/d)']       # these should cancel out
        },
'PON' : {},
'N-Diat' : {
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
'DetN' :    {'DetN,dBurDetN' : ['DetNS1,dBurS1DetN', 'DetNS2,dBurS1DetN'],
             'DetN,dDigDetN' : ['DetNS1,dDigS1DetN', 'DetNS2,dDigS1DetN']},
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

}






