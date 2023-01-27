using Profile, PProf

include("runner.jl")

function bro()
	print("Bro\n")
end

#=
function running()
	julia("./runner.jl")
end
=#

#@profile bro()  
#print(pprof())


@profile run_ABM()
print(pprof())
