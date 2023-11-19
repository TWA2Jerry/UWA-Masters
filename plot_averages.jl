using Plots
include("io_file.jl")

mdf_file = open("mdf_file.txt", "r")
adf_file = open("adf_file.txt", "r")

agent_lines = readlines(adf_file)
model_lines = readlines(mdf_file)

parsed_model_lines = 

