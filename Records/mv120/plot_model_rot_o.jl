#include("../../load_randstep.jl")
include("prog.h")
#model = initialise(pos_vels_file, 5900)
model =  AgentsIO.load_checkpoint("./Records/mv120/mv120_s5900")

fig, ax = return_thesis_figures(model, colourmap_arg = :cool, fig_box = ((4000, 3000), (8000, 7000)), colourbarlabel_arg = L"\Phi_{R*}", marker = arrow_marker, marker_size = 5)
