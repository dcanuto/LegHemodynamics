function applyperipheryics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    numupper = 0;
    if restart == "yes"
        branches = old["branches"];
    end
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            if restart == "no"
                system.branches.term[i].P[1,1] = (
                    system.branches.P[i][1,system.solverparams.JL]*mmHgToPa);
                system.branches.term[i].V[1,1] = (
                system.branches.term[i].P[1,1]*system.branches.term[i].C[1] +
                system.branches.term[i].V0[1]);
                system.branches.term[i].P[1,2] = 80*mmHgToPa;
                system.branches.term[i].V[1,2] = (
                system.branches.term[i].P[1,2]*system.branches.term[i].C[2] +
                system.branches.term[i].V0[2]);
                system.branches.term[i].P[1,3] = 40*mmHgToPa;
                system.branches.term[i].V[1,3] = (
                system.branches.term[i].P[1,3]*system.branches.term[i].C[3] +
                system.branches.term[i].V0[3]);
                system.branches.term[i].P[1,4] = 5*mmHgToPa;
                # system.branches.term[i].P[1,4] = 4.9*mmHgToPa;
                system.branches.term[i].V[1,4] = (
                system.branches.term[i].P[1,4]*system.branches.term[i].C[4] +
                system.branches.term[i].V0[4]);
                # system.branches.term[i].P[1,5] = 4.5*mmHgToPa;
                system.branches.term[i].P[1,5] = 4.5*mmHgToPa;
                system.branches.term[i].V[1,5] = (
                system.branches.term[i].P[1,5]*system.branches.term[i].C[5] +
                system.branches.term[i].V0[5]);
            elseif restart == "yes"
                term = branches["term"][i];
                system.branches.term[i].P[1,1] = term["P"][end,1];
                system.branches.term[i].V[1,1] = term["V"][end,1];
                system.branches.term[i].Q[1,1] = term["Q"][end,1];
                system.branches.term[i].P[1,2] = term["P"][end,2];
                system.branches.term[i].V[1,2] = term["V"][end,2];
                system.branches.term[i].Q[1,2] = term["Q"][end,2];
                system.branches.term[i].P[1,3] = term["P"][end,3];
                system.branches.term[i].V[1,3] = term["V"][end,3];
                system.branches.term[i].Q[1,3] = term["Q"][end,3];
                system.branches.term[i].P[1,4] = term["P"][end,4];
                system.branches.term[i].V[1,4] = term["V"][end,4];
                system.branches.term[i].Q[1,4] = term["Q"][end,4];
                system.branches.term[i].P[1,5] = term["P"][end,5];
                system.branches.term[i].V[1,5] = term["V"][end,5];
                system.branches.term[i].Q[1,5] = term["Q"][end,5];
            end
            if system.branches.group[i] == "upper"
                numupper+=1
            end
        end
    end

    if numupper > 0
        if restart == "no"
            system.svc.P[1,1] = 4.56*mmHgToPa;
            system.svc.V[1,1] = system.svc.P[1,1]*system.svc.C + system.svc.V0;
            system.ivc.P[1,1] = 4.52*mmHgToPa;
            system.ivc.V[1,1] = system.ivc.P[1,1]*system.ivc.C + system.ivc.V0;
        elseif restart == "yes"
            svc = old["svc"];
            ivc = old["ivc"];
            system.svc.P[1,1] = svc["P"][end,1];
            system.svc.V[1,1] = svc["V"][end,1];
            system.svc.Q[1,1] = svc["Q"][end,1];
            system.ivc.P[1,1] = ivc["P"][end,1];
            system.ivc.V[1,1] = ivc["V"][end,1];
            system.ivc.Q[1,1] = ivc["Q"][end,1];
        end
    else
        if restart == "no"
            system.ivc.P[1,1] = 3*mmHgToPa;
            system.ivc.V[1,1] = system.ivc.P[1,1]*system.ivc.C + system.ivc.V0;
        elseif restart == "yes"
            ivc = old["ivc"];
            system.ivc.P[1,1] = ivc["P"][end,1];
            system.ivc.V[1,1] = ivc["V"][end,1];
            system.ivc.Q[1,1] = ivc["Q"][end,1];
        end
    end
end
