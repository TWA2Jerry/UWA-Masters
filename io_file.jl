#=
###Initialise the model
model = initialise(1000*sqrt(12), 1)
print("Number of agents is $(nagents(model))\n")

figure, _ = abmplot(model)
figure # returning the figure displays it
save("shannon_flock.png", figure)
=#     

#=
compac_frac_file = open("compaction_frac.txt", "w")
mean_a_file = open("mean_area.txt", "w")
rot_o_file = open("rot_order.txt", "w")
rot_o_alt_file = open("rot_order_alt.txt", "w")
mean_speed_file = open("mean_speed.txt", "w")
=#

#=
###Animate
model = initialise(1000.0*sqrt(12), 1);
#ac(agent) =  get(balance, abs(agent.A-model.target_area)/(pi*rho^2))
plotkwargs = (; ac = get(balance, 0.7), as  = 10, am = :diamond)

print("Gone past the thang")

abmvideo(
    "Colour_Test.mp4", model, agent_step!, model_step!;
    spf = 1,
        framerate = 12, frames = 24,
    title = "Shannon flock",
        showstep = true,
        ac = get(balance, 0.7), as = 10, am = :diamond
)

print("Finished the vid\n")
=#


function do_io_stuff(compac_frac_file, mean_a_file, rot_o_file, rot_o_alt_file, mean_speed_file)
	
        close(compac_frac_file)
        close(mean_a_file)
        close(rot_o_file)
        close(rot_o_alt_file)
        close(mean_speed_file)

        compac_frac_file = open("compaction_frac.txt", "r")
        mean_a_file = open("mean_area.txt", "r")
        rot_o_file = open("rot_order.txt", "r")
        rot_o_alt_file = open("rot_order_alt.txt", "r")
        mean_speed_file = open("mean_speed.txt", "r")

cf_array = []
ma_array = []
rot_o_array = []
rot_o_alt_array = []
ms_array = []

for i in 0:no_steps
        push!(cf_array, [])
        push!(ma_array, [])
        push!(rot_o_array, [])
        push!(rot_o_alt_array, [])
        push!(ms_array, [])
end

cf_lines = readlines(compac_frac_file)
ma_lines = readlines(mean_a_file)
rot_o_lines = readlines(rot_o_file)
rot_o_alt_lines = readlines(rot_o_alt_file)
ms_lines = readlines(mean_speed_file)

print("The first thing read from the compac_frac_file was $(cf_lines[1])\n")
for line in cf_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                #print("The element read was $(split_line[i])\n")
                push!(cf_array[i], split_line[i])
        end
end

for line in ma_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(ma_array[i], split_line[i])
        end
end

for line in rot_o_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(rot_o_array[i], split_line[i])
        end
end

for line in rot_o_alt_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(rot_o_alt_array[i], split_line[i])
                #print("The split line was $(split_line[i])\n")
        end
end

for line in ms_lines
        split_line = parse.(Float64, split(line, " "))
        for i in 1:length(split_line)
                push!(ms_array[i], split_line[i])
                #print("The split line was $(split_line[i])\n")
        end
end

cf_ave_file = open("cf_ave.txt", "w")
ma_ave_file = open("ma_ave.txt", "w")
rot_o_ave_file = open("rot_o_ave.txt", "w")
rot_o_alt_ave_file = open("rot_o_alt_ave.txt", "w")
mean_speed_file = open("mean_speed_ave.txt", "w")

for i in 0:no_steps
        write(cf_ave_file, "$i $(mean(cf_array[i+1]))\n")
        write(ma_ave_file, "$i $(mean(ma_array[i+1]))\n")
        write(rot_o_ave_file, "$i $(mean(rot_o_array[i+1]))\n")
        write(rot_o_alt_ave_file, "$i $(mean(rot_o_alt_array[i+1]))\n")
        write(mean_speed_file, "$i $(mean(ms_array[i+1]))\n")
end

end
