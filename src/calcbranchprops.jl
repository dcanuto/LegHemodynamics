function calcbranchprops!(system::LegSystem,old=Dict("a"=>0),restart="no")
    for i in 1:length(system.branches.ID)
        # Moens-Korteweg wave speed, m/s
        push!(system.branches.c0,[sqrt(0.5*system.branches.beta[i][end]/
            system.solverparams.rho)*system.branches.A0[i]^0.25])
    end
    if restart == "no"
        push!(system.branches.termscalings,1.)
        push!(system.branches.termscalings,1.)
        push!(system.branches.termscalings,1.)
        push!(system.branches.termscalings,1.)
        push!(system.branches.termscalings,1.)
    elseif restart == "yes"
        append!(system.branches.termscalings,old["termscalings"])
    end
end
