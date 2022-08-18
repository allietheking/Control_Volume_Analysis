

# To create the balance tables, you can either run the scripts one at a time, or you can use this script to run them all at once.
# To specify the run you want to process, the directory where you want the balance table output, and other user input variables, 
# make changes to the step0_config.py script. The configuration variables specified in step0_config.py are imported by each of the 
# following modules
#
# created by Allie August 2022



python step1_create_balance_tables
python step2_create_balance_tables_for_composite_parameters
python step3_group_reactions_for_composite_parameters
python step4_check_mass_conservation
python step5_compile_balance_tables_into_groups
python step6_aggregate_in_time