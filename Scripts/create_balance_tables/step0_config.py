
''' Configuration module for ALL of the balance table scripts 
This is where we set global variables that are used by all the scripts. 
Many other settings are specified in the individual scripts as well, in the 
USER INPUT section, so make sure to take a look there as well to see what the
scripts are doing

alliek august 2022
'''

# this is the run you want to process
runid = 'G141_13to18_207'

# this is the root directory for model input and output
# the scripts look in this directory for the shapefiles used to specify the run, the run folder containing the
model_inout_dir = '/richmondvol1/hpcshared/'

# is this a delta run? (note that if so, scripts assume it is full resolution)
is_delta = False

# base level substances to process (either 'all' or a list of substances -- warning, processing all of them takes a long time)
# substance_list = 'all'
substance_list = ['continuity', 'don', 'detns1', 'detns2', 'diat', 'diats1', 'green', 'nh4', 'no3', 'oxy', 'zoopl_e', 'zoopl_r', 'zoopl_v']

# list of time averaging schemes to apply (saves space to skip some if we don't need them)
tavg_list = ['Cumulative', 'Filtered', 'Annual', 'Seasonal', 'Monthly', 'Weekly']

# directory for saving the balance tables (could be in main repo or in run folder)
if is_delta:
	balance_table_dir = '/chicagovol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/Delta_%s' % runid
else:
	balance_table_dir = '/chicagovol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Balance_Tables/%s' % runid

# directory where aggregated groups are defined (this is an input directory, probably don't want to change it)
group_def_dir = '/chicagovol1/hpcshared/NMS_Projects/Control_Volume_Analysis/Definitions'

# path to stompy directory (used in step 1)
stompy_dir = '/opt/software/rusty/stompy/newest_commit/stompy'

# float format for csv files
float_format = '%1.16e'

# set tolerance for mass conservation error as a percentage of whatever variable we chose to normalize by
error_tol_percent = 0.1

# abort if error tolerance is exceeded???
abort_for_mass_cons_error = False


