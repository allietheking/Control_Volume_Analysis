

#################################################################
# IMPORT PACKAGES
#################################################################

import sys,os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt



#################################################################
# USER INPUT
#################################################################

runid_list = ['FR13_021','FR17_014','FR18_004','FR13_003','FR17_003', 'G141_13to18_125']
group_list = ['Whole_Bay','A']

sort_rx_names = True
sort_order = ['NH4','NO3','PON','N-D','N-A','THI','EAC','SHO']

composite_parameter = 'TN'

balance_table_dir = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Balance_Tables' # balance tables should be organized by runid within this folder

plot_dir = r'.\plots'

color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',
              u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',
              u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

style_list = ['-','-','-','-','-','-','-','-','-','-',
              '--','--','--','--','--','--','--','--','--','--',
              ':',':',':',':',':',':',':',':',':',':']

#################################################################
# MAIN
#################################################################

  
# run id
for runid in runid_list: 

    # group
    for group in group_list:
        
        # input file -- balance table for raw parameters, aggregated into groups
        input_balance_table_filename = os.path.join(balance_table_dir,runid,'Balance_Table_By_Group_Composite_Parameter_%s.csv' % composite_parameter)
        
        # output plot file name 
        output_plot_filename = os.path.join(plot_dir, '%s_%s_%s_Rx.png' % (runid, group, composite_parameter))
        
        # read balance table
        df = pd.read_csv(input_balance_table_filename) 
                
        # select group
        ind = df['group'] == group
        df = df.loc[ind]
        
        # compute time
        time = df['time'].astype('datetime64[ns]').values
        
        # pick out the net reactions
        netrx = df['%s,Net Reaction (Mg/d)' % composite_parameter].values

        # get names of reactions
        rxnames = list(df.columns[15:])
        
        # sort the reaciton terms according to first three characters
        if sort_rx_names:
            rxnames_sorted = []
            rxnames_unsorted = rxnames.copy()
            for sort_item in sort_order:
               for rx in rxnames:
                    if sort_item == rx[0:3]:
                        rxnames_sorted.append(rx)
                        rxnames_unsorted.remove(rx)
            for rx in rxnames_unsorted:
               rxnames_sorted.append(rx)
            rxnames = rxnames_sorted.copy()


        
        # plot
        fig, ax = plt.subplots(figsize=(18,10))
        ax.plot(time,netrx,label='NET REACTION (Mg/d)', color='k', linewidth=4)
        for i in range(len(rxnames)):
            ax.plot(time,df[rxnames[i]],label=rxnames[i],color=color_list[i], linestyle=style_list[i], linewidth=2)
        ax.legend()
        ax.set_title('Runid = %s, Composite Parameter = %s, Group = %s' % (runid,composite_parameter,group), fontsize=16)
        fig.tight_layout()
        fig.savefig(output_plot_filename)
        plt.close('all')