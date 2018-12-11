function calcbranchprops!(system::LegSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        for i in 1:length(system.branches.ID)
            # Moens-Korteweg wave speed, m/s
            push!(system.branches.c0,[sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*system.branches.A0[i]^0.25])
        end
    elseif restart == "yes"
        for i = 1:length(system.branches.ID)
            push!(system.branches.c0,[old["c0"][i][end]]);
        end
    end
end
