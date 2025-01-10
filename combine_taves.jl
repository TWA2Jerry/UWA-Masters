##Let mdf be the total data frame for all steps, all simulations, all target areas. Here we sort that data frame into groups with a group for each tdod and simulation number. Then we combine the groups, averaging across time steps. The resulting df has average rot o for each sim and tdod. Then, we group by tdods again, and combine again, thus taking the average across sims.  
tdod_grouped_mdf = groupby(mdf, [:target_area_arg, :simulation_number_arg])
simulation_combined_mdf = combine(tdod_grouped_mdf, :target_area_arg => mean, :rot_o_alt => mean)
tdod_simulation_grouped_mdf = groupby(simulation_combined_mdf, :target_area_arg)
tdod_combined_mdf = combine(tdod_simulation_grouped_mdf, :rot_o_alt_mean => mean, :rot_o_alt_mean => std)


