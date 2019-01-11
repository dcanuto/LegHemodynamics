function applyperipheryics!(system::LegSystem,old=Dict("a"=>0),restart="no")
    numupper = 0;
    if restart == "yes"
        branches = old["branches"];
    end
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            if restart == "no"
                system.branches.term[i].P[1,1] = (
                    system.branches.P[i][system.solverparams.JL,1]*mmHgToPa);
            elseif restart == "yes"
                term = branches["term"][i];
                system.branches.term[i].P[1,1] = term["P"][end,1];
                system.branches.term[i].Q[1,1] = term["Q"][end,1];
            end
        end
    end
end
