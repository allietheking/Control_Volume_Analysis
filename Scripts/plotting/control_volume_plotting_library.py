#/chicagovol1/hpcshared/open_bay/bgc/figures


#/richmondvol1/hpcshared/Grid141/WY13to18/G141_13to18_140
#/richmondvol1/hpcshared/Grid141/WY2013/G141_13_016
#/richmondvol1/hpcshared/Full_res/WY2013/FR13_025


import numpy as np

def make_concise_runid_list_string(runid_list):

	'''Given a list of run IDs for an NMS open bay biogeochemical model run, 
	return a concise string listing all the run IDs, intended for naming postprocessing
	figures and the directories where we will store them.

	Usage:

		concise_runid_list_string = make_concise_runid_list_string(runid_list)

	An example valid runid_list is 

		['FR18_006', 'G141_17_017','G141_13to18_207','G141_13to18_197','G141_13_016', 'G141_13_017','FR13_025']

	and the value of concise_runid_list_string returned for this runid_list would be

	 	'AGG13to18_207_197_AGG13_016_017_AGG17_017_FR13_025_FR18_006'
	
	'''

	# replace 'G141_' with 'AGG' in runid's
	for i,runid in enumerate(runid_list):
		runid_list[i] = runid.replace('G141_','AGG')
	
	# come up with a corresponding list of just the run prefixes and just the run numbers
	run_pref_list = []
	run_numb_list = []
	for runid in runid_list:
		run_pref_list.append(runid.split('_')[0])
		run_numb_list.append(runid.split('_')[1])
	
	# come up with a sorted list of unique run prefixes (such as AGG13to18, AGG17, FR17) for the AGG and FR runs separately
	run_pref_unique_AGG = []
	run_pref_unique_FR = []
	for runid in runid_list:
		if 'AGG' in runid:
			run_pref_unique_AGG.append(runid.split('_')[0])
		else:
			run_pref_unique_FR.append(runid.split('_')[0])
	run_pref_unique_AGG = np.sort(np.unique(run_pref_unique_AGG))
	run_pref_unique_FR = np.sort(np.unique(run_pref_unique_FR))
	
	# move any runs with 'to' in the runid up to the front
	def move_to_up(run_pref_unique):
		i_move = []
		i_stay = []
		for i in range(len(run_pref_unique)):
			if 'to' in run_pref_unique[i]:
				i_move.append(i)
			else:
				i_stay.append(i)
		run_pref_unique = np.concatenate((run_pref_unique[i_move], run_pref_unique[i_stay]))
		return run_pref_unique
	run_pref_unique_AGG = move_to_up(run_pref_unique_AGG)
	run_pref_unique_FR = move_to_up(run_pref_unique_FR)
	
	# concatenate FR and AGG runs to make one list of unique runid prefixes
	run_pref_unique = np.concatenate((run_pref_unique_AGG, run_pref_unique_FR))
	
	# now loop through the sorted unique run prefixes, find all the runs with that prefix, get the corresponding run numbers,
	# sort the run numbers, and append to the folder name
	concise_runid_list_string = ''
	for run_pref in run_pref_unique:
		concise_runid_list_string = concise_runid_list_string + run_pref + '_'
		i_list = np.where(run_pref==np.array(run_pref_list))[0]
		for i in i_list:
			concise_runid_list_string = concise_runid_list_string + run_numb_list[i] + '_'
	concise_runid_list_string = concise_runid_list_string[0:-1]
	
	return concise_runid_list_string