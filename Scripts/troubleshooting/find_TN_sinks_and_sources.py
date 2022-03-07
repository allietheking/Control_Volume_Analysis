

#################################################################
# IMPORT PACKAGES
#################################################################

import sys
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import datetime as dt
import itertools


#################################################################
# USER INPUT
#################################################################
  
# composite parameters -- for each parameter, include a stoiciometric multiplier;
# for example for diatoms, the fraction contributing to nitrogen is 0.14,
# so when computing total nitrogen, need to multiply diatom masses by 0.14
composite_parameters = {'TN'  : [[1,'NH4'], [1,'NO3'], [1,'DON'], [1,'PON1'], [0.15,'Diat'], [0.1818,'Zoopl_V'], [0.1818,'Zoopl_R'], [0.1818,'Zoopl_E']]}
 
# input file -- balance table for raw parameters, aggregated into groups
input_balance_table_filename = 'Balance_Tables/Balance_Table_All_Parameters_By_Group.csv'

# define zero
zero_cut = 0.00001

#########################
# MAIN
#########################

# read main balance table with raw parameters and select only group A
df_master = pd.read_csv(input_balance_table_filename)
df_master = df_master.loc[df_master['group']=='A']

# get list of constituent parameters and corresponding multipliers
mult_list = []
param_list = []
for mp in composite_parameters['TN']:
    mult_list.append(mp[0])
    param_list.append(mp[1])

# make a list of columns that includes the groups, time stamps, and all the reaction
# terms corresponding to the constituent parameters
columns = []
for col in df_master.columns:
    if ',d' in col:
        if not ',dMass/dt' in col:
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

## delete the columns we already know to be sources and sinks
#df.drop(columns=['NH4,dMinDetNS1 (Mg/d)',
#                 'PON1,dResS1DetN (Mg/d)',
#                 'PON1,dSedPON1 (Mg/d)',
#                 'NO3,dDenitWat (Mg/d)',
#                 'NO3,dDenitSed (Mg/d)'], inplace=True)
# delete the columns we already know to be sources and sinks
df.drop(columns=['NH4,dMinDetNS1 (Mg/d)', # source for TN
                 'NH4,dMinDetNS2 (Mg/d)', # source for TN
                 'NH4,dNH4Upt (Mg/d)', # part of uptake balance
                 'NO3,dDenitWat (Mg/d)', # sink for TN
                 'NO3,dDenitSed (Mg/d)', # sink for TN
                 'NO3,dNO3Upt (Mg/d)', # part of uptake balance
                 'PON1,dSedPON1 (Mg/d)', # sink for TN 
                 'PON1,dResS1DetN (Mg/d)', # source for TN
                 'PON1,dResS2DetN (Mg/d)', # source for TN
                 'Diat,dSedDiat (Mg/d)', # flagged by exploratory algorithm v1
                 'Diat,dPPDiat (Mg/d)', # part of uptake balance
                 'Diat,dcPPDiat (Mg/d)'], # part of uptake balance
                 inplace=True) 
columns = list(df.columns)
Nc = len(columns)

# check if the remaining columns sum to zero
if not any(df.sum(axis=1).abs()>zero_cut):
    print('residual columns sum to zero -- do not need to seek more sources/sinks')

# try eliminating sets of n reaction terms for n = 1, 2, 3, ... 
n = -1
while n < 6:
    n = n + 1
    print('Trying sets of additional sources/sinks of length n = %d' % n)
    for col_ind_set in itertools.combinations(range(Nc), n):
        df1 = df.copy(deep=True)
        dropcols = []
        for col_ind in col_ind_set:
            dropcols.append(columns[col_ind])
        df1.drop(columns=dropcols, inplace=True)
        if not any(df1.sum(axis=1).abs()>zero_cut):
            print('Here are the missing sources/sinks:')
            print(dropcols)
            sys.exit()
            
            
        
        
    
