function discretizeperiphery!(system::CVSystem)
    for i in 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            system.branches.term[i].P = zeros(system.solverparams.numsteps+1,5);
            system.branches.term[i].V = zeros(system.solverparams.numsteps+1,5);
            system.branches.term[i].Q = zeros(system.solverparams.numsteps+1,5);
        end
    end

    append!(system.svc.P,zeros(system.solverparams.numsteps+1))
    append!(system.svc.V,zeros(system.solverparams.numsteps+1))
    append!(system.svc.Q,zeros(system.solverparams.numsteps+1))
    append!(system.ivc.P,zeros(system.solverparams.numsteps+1))
    append!(system.ivc.V,zeros(system.solverparams.numsteps+1))
    append!(system.ivc.Q,zeros(system.solverparams.numsteps+1))
end
