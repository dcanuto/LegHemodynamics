function applybranchics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        Pstart = 70*mmHgToPa;
        for i in 1:length(system.branches.ID)
            system.branches.A[i][:,1] = (Pstart/system.branches.beta[i][end] +
                sqrt(system.branches.A0[i][end]))^2;
            system.branches.Q[i][:,1] = 0;
            system.branches.P[i][:,1] = Pstart/mmHgToPa;
        end
    elseif restart == "yes"
        branches = old["branches"];
        for i in 1:length(system.branches.ID)
            if system.solverparams.JL == old["solverparams"]["JL"]
                system.branches.A[i][:,1] = branches["A"][i][:,end];
                system.branches.Q[i][:,1] = branches["Q"][i][:,end];
                system.branches.P[i][:,1] = branches["P"][i][:,end];
            else
                A = branches["A"][i][:,end];
                Q = branches["Q"][i][:,end];
                P = branches["P"][i][:,end];
                Aitp = Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
                Qitp = Interpolations.interpolate(Q, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
                Pitp = Interpolations.interpolate(P, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
                x = 0:1:(old["solverparams"]["JL"]-1);
                x = x*branches["k"][1];
                sAitp = scale(Aitp,x);
                sQitp = scale(Qitp,x);
                sPitp = scale(Pitp,x);
                xq = 0:1:(system.solverparams.JL-1);
                xq = xq*system.branches.k[i];
                Aq = [sAitp[j] for j in xq];
                Qq = [sQitp[j] for j in xq];
                Pq = [sPitp[j] for j in xq];
                system.branches.A[i][:,1] = Aq;
                system.branches.Q[i][:,1] = Qq;
                system.branches.P[i][:,1] = Pq;
            end
        end
    end
end
