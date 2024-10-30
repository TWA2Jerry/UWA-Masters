function return_stability_rot_o(data)
	filtered_data = filtered_data = filter(row -> row.step >= 55000, mdf)
	stability_rotation = calculate_stability_rotation(filtered_data)
	single_cluster_rot_o = calculate_single_cluster_rot_o(filtered_data)
	stability_score = calculate_stability_score(filtered_data)
	return stability_rotation, single_cluster_rot_o, stability_score
end
