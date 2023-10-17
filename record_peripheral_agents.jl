periphery_file = open("periphery_file.txt", "w")

function detect_periphery(cell, model_n)
        on_periphery::Int32 = 0
        for vertex in cell
                if(vertex[2] == 0 || vertex[3] == 0)
                        on_periphery = 1
                        break
                end
        end
	return on_periphery
end


function detect_write_periphery(area, cell, model_n)
	on_periphery::Int32 = 0
	for vertex in cell
		if(vertex[2] == 0 || vertex[3] == 0)
			on_periphery = 1
			break
		end
	end
	write(periphery_file, "$area $on_periphery\n")
end


