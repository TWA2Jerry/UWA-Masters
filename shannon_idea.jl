//Add no_positions to agent definition parameters

num_positions::Int32 = 0

//Then throughout your move gradient code, if it moves past the valid stage, and into the calculation
num_positions += 1

//Then after all the thingos have been done
agent.num_positions = num_positions
