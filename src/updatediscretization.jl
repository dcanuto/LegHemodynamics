function updatediscretization!(system::LegSystem)
    # determine time shift
    system.solverparams.tshift = system.t[end]-
        sum(system.solverparams.th*system.solverparams.numbeats);

    # update total number of time steps
    ntoadd = ceil((system.solverparams.th-
        system.solverparams.tshift)/system.solverparams.h);
    system.solverparams.numsteps+=ntoadd;

    # update discrete times
    ttoadd = [0:1:ntoadd;]*system.solverparams.h + system.t[end];
    append!(system.t,ttoadd[2:end]);

    # allocate additional space for solution variables
    for i = 1:length(system.branches.ID)
        system.branches.A[i] = [system.branches.A[i] zeros(system.solverparams.JL,length(ttoadd[2:end]))];
        system.branches.Q[i] = [system.branches.Q[i] zeros(system.solverparams.JL,length(ttoadd[2:end]))];
        system.branches.P[i] = [system.branches.P[i] zeros(system.solverparams.JL,length(ttoadd[2:end]))];
        if isempty(system.branches.children[i])
            system.branches.term[i].P = [system.branches.term[i].P;zeros(length(ttoadd[2:end]),1)];
            system.branches.term[i].Q = [system.branches.term[i].Q;zeros(length(ttoadd[2:end]),1)];
        end
    end

end
