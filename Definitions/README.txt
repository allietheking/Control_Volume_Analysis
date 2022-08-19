
Here's a description of everything in the various subfolders in the Definitions folder (same directory as 
this README) Here you can find shapefiles and figures illustrating the regions represented by the polygons and 
transects you'll find in the BGC model *.his and *-his.bal output files. You'll also find the definitions
and illustrations of larger "groups" comprised of multiple polygons, such as the whole bay, or the south 
bay (RMP definition). Some of these groups were the brainchild of Dave, and are illustrated in this 
presentation: https://docs.google.com/presentation/d/178gRZfdF5x9YbxwKmeAvr-Dsi2mk_ug6UK2PKYQxNPw/edit#slide=id.p

-----------------------
model_input_shapefiles
-----------------------

In our BGC model setup scripts, we read in shapefiles to define "monitoring areas" and "transects".
We have saved a copy of these shapefiles here in the "model_input_shapefiles" folder, and in this folder
you will also find *.png files with the same name as the shapefiles. Here's a guide:

Agg_mod_contiguous.shp defines the "monitoring areas" for the full resolution run. This *.his file gives
concentrations averaged across these regions, and the *-his.bal file gives the rates of all the reactions as
well as the net loading and the net transport into these regions. Stompy gives the monitoring areas the 
the names "polygon0", "polygon1", "polygon2", etc. in the *.his and *-bal.his file. The number corresponding 
to each polygon is illustrated in figure Agg_mod_contiguous.png, so take a look at that figure to see for 
example what "polygon23" corresponds to. 

Agg_exchange_lines.shp defines the "transects" between the monitoring areas for the full resolution runs. The 
*.his file reports the net flux across each transect during the *.his time step. Stompy gives the transects 
the names "transect0", "transect1", "transect2", etc. and the figure Agg_exchange_lines.png illustrates the 
transect corresponding to each number. 

Similarly, for the aggregated grid model, the shapefiles Agg_mod_contiguous_141.shp and Agg_exchange_lines_141.shp
define the monitoring area and transects, and the *.png figures of the same name illustrate where each polygon
and transect number are located.

-----------------------
group_definitions
-----------------------

The model output is at the level of the polygons and transects defined in "model_input_shapefiles". 
We wanted to define larger regions, such as the whole bay, or different subembayments, or the west shoal
of south bay, comprised of multiple polygons, and we wanted to compute the flux across the N, S, E, and 
W sides of these regions. In the group_definitinos folder you will find definintions of these larger groups 
and their N, S, E, and W faces for both full resoluton and aggregateed grid runs.

-----------------------
group_shapefiles
-----------------------

This folder contains shapefiles for the "groups" and their connectivites. You will also find in this folder the
script that created these shapefiles from the group definitions found in the "group_definitions" folder. Some of the
plotting scripts rely on these shapefiles

-----------------------
plots_of_groups
-----------------------

The plots in this folder illustrate the groups and their connectivity as defined in the "group_definitions" folder

-----------------------
nice_group_plots
-----------------------

These plots of various sets of groups, such as the subembayments, or the south bay, are intended for 
presentations and manuscripts
