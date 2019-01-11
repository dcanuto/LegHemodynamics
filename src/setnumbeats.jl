function setnumbeats!(system::LegSystem,n::Int64)
    # total time needed for all past/current cardiac cycles
    ttotal = system.solverparams.th*(system.solverparams.numbeats+1);

    # increment number of cycles if new cycle is starting
    if system.t[n+2] >= ttotal
        println("Old number of completed cardiac cycles: $(system.solverparams.numbeats)")
        println("New number of completed cardiac cycles: $(system.solverparams.numbeats+1)")
        println("Total number of cycles to complete: $(system.solverparams.numbeatstotal)")
        system.solverparams.numbeats+=1;
    end
end
