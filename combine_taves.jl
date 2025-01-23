##Let mdf be the total data frame for all steps, all simulations, all target areas. Here we sort that data frame into groups with a group for each tdod and simulation number. Then we combine the groups, averaging across time steps. The resulting df has average rot o for each sim and tdod. Then, we group by tdods again, and combine again, thus taking the average across sims.  
using CSV
using DataFrames
using StatsBase
record = "stdod73" 
mdf = DataFrame(CSV.File("Records/$(record)/$(record)_1stdod_mdata"))
if(columnindex(mdf, :q_arg) > 0)
	mdf = select!(mdf, Not([:q_arg]))
end

if(columnindex(mdf, :m_arg) > 0)
    mdf = select!(mdf, Not([:m_arg]))
end 

if(columnindex(mdf, :fov_arg) > 0)
    mdf = select!(mdf, Not([:fov_arg]))
end 


for i in [10, 100, 1000, 2000]
	mdft = DataFrame(CSV.File("Records/$(record)/$(record)_$(i)stdod_mdata"))	
	if(columnindex(mdft, :q_arg) > 0)
		mdft = select!(mdft, Not([:q_arg]))
	end 

	if(columnindex(mdft, :m_arg) > 0)	
		mdft = select!(mdft, Not([:m_arg]))
	end

	if(columnindex(mdft, :fov_arg) > 0)
		mdft = select!(mdft, Not([:fov_arg]))
	end
	append!(mdf, mdft)
end

mdft = DataFrame(CSV.File("Records/$(record)/$(record)_$(22000)_mdata"))
if(columnindex(mdft, :q_arg) > 0)
    mdft = select!(mdft, Not([:q_arg]))
end

if(columnindex(mdft, :m_arg) > 0)
    mdft = select!(mdft, Not([:m_arg]))
end

if(columnindex(mdft, :fov_arg) > 0)
    mdft = select!(mdft, Not([:fov_arg]))
end

append!(mdf, mdft)


tdod_grouped_mdf = groupby(mdf, [:target_area_arg, :simulation_number_arg])
simulation_combined_mdf = combine(tdod_grouped_mdf,  :target_area_arg => mean, :mean_speed => mean, :rot_o_alt => mean)
tdod_simulation_grouped_mdf = groupby(simulation_combined_mdf, :target_area_arg)
tdod_combined_mdf = combine(tdod_simulation_grouped_mdf, :mean_speed_mean => mean, :mean_speed_mean => std, :rot_o_alt_mean => mean, :rot_o_alt_mean => std)


