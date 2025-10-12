#use this in julia interactive mode while in the UWA Masters folder. 
include("load_randstep.jl")
pos_vels_file = open("mv133pos_vels.txt", "r")
model = initialise(pos_vels_file, 1713)
fig, ax = return_thesis_figures(model, fig_box = ((1500, 0), (2500, 1000)), marker = arrow_marker, marker_size= 10, colourmap_arg = :cool, colourbarvisible_arg = 0)
