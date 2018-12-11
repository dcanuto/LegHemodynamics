function calcbranchprops!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        for i in 1:length(system.branches.ID)
            # reference area, m^2
            push!(system.branches.A0,[π*(system.branches.radiusincm[i]*cmTom)^2])
            # stiffness parameter, Pa/m
            push!(system.branches.beta,[(sqrt(π)*system.branches.thicknessincm[i]*
                cmTom*system.branches.YoungsModinMPa[i]*MPaToPa)/
                ((1-system.solverparams.nu^2)*system.branches.A0[i][end])])
            # Moens-Korteweg wave speed, m/s
            push!(system.branches.c0,[sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*system.branches.A0[i][end]^0.25])
        end
    elseif restart == "yes"
        for i = 1:length(system.branches.ID)
            push!(system.branches.A0,[old["A0"][i][end]]);
            push!(system.branches.beta,[old["beta"][i][end]]);
            push!(system.branches.c0,[old["c0"][i][end]]);
        end
    end
end
