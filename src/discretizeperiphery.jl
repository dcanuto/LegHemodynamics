function discretizeperiphery!(system::LegSystem)
    for i in 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            system.branches.term[i].P = zeros(system.solverparams.numsteps+1,1);
            system.branches.term[i].Q = zeros(system.solverparams.numsteps+1,1);
        end
    end
end
