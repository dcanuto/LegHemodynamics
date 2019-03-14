function rediscretizet!(system::LegSystem,n::Int64)
    # time remaining
    trem = system.solverparams.th*(system.solverparams.numbeats+1) - system.t[n+1];
    # number of steps remaining
    nrem = ceil(trem/system.solverparams.h);
    # update numsteps
    system.solverparams.numsteps = n+nrem;
    # update discrete times
    ttoadd = [0:1:nrem;]*system.solverparams.h + system.t[n+1];
    resize!(system.t,n+1)
    append!(system.t,ttoadd[2:end]);
    # resize allocated space for solution variables
    for i = 1:length(system.branches.ID)
        system.branches.A[i] = system.branches.A[i][:,1:n+1];
        system.branches.A[i] = [system.branches.A[i] zeros(system.solverparams.JL,nrem)];
        system.branches.Q[i] = system.branches.Q[i][:,1:n+1];
        system.branches.Q[i] = [system.branches.Q[i] zeros(system.solverparams.JL,nrem)];
        system.branches.P[i] = system.branches.P[i][:,1:n+1];
        system.branches.P[i] = [system.branches.P[i] zeros(system.solverparams.JL,nrem)];
        if isempty(system.branches.children[i])
            system.branches.term[i].P = system.branches.term[i].P[1:n+1];
            system.branches.term[i].P = [system.branches.term[i].P;zeros(length(ttoadd[2:end]),1)];
            system.branches.term[i].Q = system.branches.term[i].Q[1:n+1];
            system.branches.term[i].Q = [system.branches.term[i].Q;zeros(length(ttoadd[2:end]),1)];
        end
    end
    return system
end
