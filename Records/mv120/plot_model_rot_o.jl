include("../../load_randstep.jl")

model = initialise(pos_vels_file, 5900)

#fig, ax = return_thesis_figures(model, colourmap_arg = :cool, colourbarlabel_arg = L"\Phi_{R*}")
fig, ax = return_thesis_figures(model, colourmap_arg = :cool, fig_box = ((2000, 3000), (7000, 8000)), colourbarlabel_arg = L"\Phi_{R*}")
