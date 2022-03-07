'''

Use new "Whole_Bay" group to work up budget

This creates a four panel plot with the following panels
- bay-wide DIN Loading
- Delta influx (net flux in from the W)
- assimilation (net Rx term)
- outflux throguh Golden Gate (net flux in from the E)

Multiple water years are plotted together, based on full resoluation runs

'''

########################################################################################
## import python packages
########################################################################################

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import datetime as dt
import matplotlib.dates as mdates
from scipy import signal
#plt.switch_backend("Agg")


#########################################################################################
## user input
#########################################################################################

wy_list = [2013,2017,2018]
run2plot_list = ['FR13_003','FR17_003','FR18_003']

## composite parameter (must match suffix of balance table)
param = 'DIN'

## output plot directory
out_dir = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Plots\Compare_Water_Years'

## here's a colorblind friendly color cycle
cfc = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# set y limits for plots
ylim_day = 700
ylim_tif = 500
ylim_cum = 55000/1000

# figure size 
fs = (8,8)

#########################################################################################
## functions
#########################################################################################

def spring_neap_filter(y, dt_days = 1.0, fcut = 1/36., N = 6):

    """
    Performs a spring-neap filter on the y series using a 6th order low pass
    butterworth filter and constant padding
   
    Usage:
       
        yf = spring_neap_filter(y, dt_days=1.0, fcut = 1/36.)
       
    Input:
       
        y    = time series
        dt_days   = time step (unit is days)
        fcut = cutoff frequency in 1/days
               
    Output:
   
        yf   = tidally filtered signal
   
    """
   
    # compute sampling frequency from time step
    fs = 1/dt_days
   
    # compute dimensionless cutoff frequency w.r.t. Nyquist frequency
    Wn = fcut/(0.5*fs)
    # compute filter coefficients
    b, a = signal.butter(N, Wn, 'low')
   
    # filter the signal
    yf = signal.filtfilt(b,a,y,padtype='constant')
   
    # return fitlered signal
    return yf

#########################################################################################
## main
#########################################################################################

# initialize figure
fig0, ax0 = plt.subplots(4,1,figsize=fs)
fig1, ax1 = plt.subplots(4,1,figsize=fs)
fig2, ax2 = plt.subplots(4,1,figsize=fs)

# dummy start time
ts = np.datetime64('2012-10-01')

for iwy in range(len(wy_list)):

    wy = wy_list[iwy]
    run2plot = run2plot_list[iwy]

    ## balance table folder
    table_dir = r'X:\hpcshared\NMS_Projects\Control_Volume_Analysis\Balance_Tables\%s' % run2plot

    ## pick the time window for cumulative plot
    t_window = np.array(['%d-10-01' % (wy-1),'%d-10-01' % wy]).astype('datetime64')

    # water year string
    water_year = 'WY%d' % wy

    # input is the composite parameter balance table, in the corresponding run folder
    input_fn = os.path.join(table_dir,'Balance_Table_By_Group_Composite_Parameter_%s.csv' % param)

    # load up the balance table data
    data = pd.read_csv(input_fn)

    # convert times from string to datetime64
    data['time'] = pd.to_datetime(data.time)

    # select the whole bay control volume
    data = data.loc[data['group']=='Whole_Bay']

    # if this is WY2018, replace the first loading value (which is zero) with the second loading value
    data['%s,Net Load (Mg/d)' % param].iloc[0] = data['%s,Net Load (Mg/d)' % param].iloc[1]

    # get subset of column names we want to apply tidal filter and cumulative sum on
    cols = data.columns[4:]

    # apply tidal filter before selecting time window to get better end effects
    data_tfi = data.copy(deep=True)
    for col in cols:
        data_tfi[col] = spring_neap_filter(data[col])

    # now pare down to time window of interest
    indt = np.logical_and(data['time']>=t_window[0], data['time']<t_window[1])
    data = data.loc[indt]
    data_tfi = data_tfi.loc[indt]

    # compute cumulutaive sum (compute time step in days)
    deltat = (data['time'].iloc[1] - data['time'].iloc[0])/np.timedelta64(1,'h')/24
    data_cum = data.copy(deep=True)
    for col in cols:
        data_cum[col] = np.cumsum(data[col]) * deltat

    # plot cumulative
    ax0[0].plot(data_cum['time']-t_window[0]+ts,data_cum['%s,Net Load (Mg/d)' % param]/1000,label=water_year)
    ax0[1].plot(data_cum['time']-t_window[0]+ts,data_cum['%s,Flux In from E (Mg/d)' % param]/1000,label=water_year)
    ax0[2].plot(data_cum['time']-t_window[0]+ts,data_cum['%s,Net Reaction (Mg/d)' % param]/1000,label=water_year)
    ax0[3].plot(data_cum['time']-t_window[0]+ts,data_cum['%s,Flux In from W (Mg/d)' % param]/1000,label=water_year)
    
    # plot tidally filtered
    ax1[0].plot(data_tfi['time']-t_window[0]+ts,data_tfi['%s,Net Load (Mg/d)' % param],label=water_year)
    ax1[1].plot(data_tfi['time']-t_window[0]+ts,data_tfi['%s,Flux In from E (Mg/d)' % param],label=water_year)
    ax1[2].plot(data_tfi['time']-t_window[0]+ts,data_tfi['%s,Net Reaction (Mg/d)' % param],label=water_year)
    ax1[3].plot(data_tfi['time']-t_window[0]+ts,data_tfi['%s,Flux In from W (Mg/d)' % param],label=water_year)

    # plot daily
    ax2[0].plot(data['time']-t_window[0]+ts,data['%s,Net Load (Mg/d)' % param],label=water_year)
    ax2[1].plot(data['time']-t_window[0]+ts,data['%s,Flux In from E (Mg/d)' % param],label=water_year)
    ax2[2].plot(data['time']-t_window[0]+ts,data['%s,Net Reaction (Mg/d)' % param],label=water_year)
    ax2[3].plot(data['time']-t_window[0]+ts,data['%s,Flux In from W (Mg/d)' % param],label=water_year)

# clean up the axes, add labels, legend, titles, grid, etc.
for (ax, units) in zip([ax0,ax1,ax2], ['\n(10$^9$g/d)', '\n(Mg/d)', '\n(Mg/d)']):
    ax[0].set_ylabel('POTW Loading' + units)
    ax[1].set_ylabel('Delta Influx' + units)
    ax[2].set_ylabel('Assimilation' + units)
    ax[3].set_ylabel('Golden Gate Outflux' + units)
    for i in range(4):
        ax[i].grid(b=True, alpha = 0.3)
        ax[i].set_xlim((ts,ts+np.timedelta64(365,'D')))
        ax[i].legend()
        ax[i].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
        ax[i].xaxis.set_major_formatter(mdates.DateFormatter('1-%b'))
ax0[0].set_title('Whole Bay Cumulative %s Budget' % param)
ax1[0].set_title('Whole Bay Tidally Fitlered %s Budget' % param)
ax2[0].set_title('Whole Bay Daily %s Budget' % param)

# create suffix with all runs in compairison
run_suff = ''
for run2plot in run2plot_list:
    run_suff = run_suff + '_' + run2plot

# format and save figures
for fig in [fig0,fig1,fig2]:
    fig.autofmt_xdate()
    fig.tight_layout()
fig0.savefig(os.path.join(out_dir, 'Whole_Bay_Budget_Cumulative_%s_%s.png' % (param, run_suff)))
fig1.savefig(os.path.join(out_dir, 'Whole_Bay_Budget_Filtered_%s_%s.png' % (param, run_suff)))
fig2.savefig(os.path.join(out_dir, 'Whole_Bay_Budget_Daily_%s_%s.png' % (param, run_suff)))


plt.close('all')