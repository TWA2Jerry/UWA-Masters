tdod_grouped_mdf = groupby(mdf, [:target_area_arg, :simulation_number_arg])
simulation_combined_mdf = combine(tdod_grouped_mdf, :target_area_arg => mean, :rot_o_alt => mean)
tdod_simulation_grouped_mdf = groupby(simulation_combined_mdf, :target_area_arg)
tdod_combined_mdf = combine(tdod_simulation_grouped_mdf, :rot_o_alt_mean => mean, :rot_o_alt_mean => std)


