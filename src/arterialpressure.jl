function arterialpressure!(system::LegSystem,n::Int64)
    for i = 1:length(system.branches.ID)
        # for j = 1:system.solverparams.JL
        #     if system.branches.A[i][j,n+2] < 0
        #         error("Area negative in pressure calculation for artery $i, node $j,
        #             time step $n. Area = $(system.branches.A[i][j,n+2]) m2.")
        #     end
        # end
        system.branches.P[i][:,n+2] .= system.branches.beta[i][end].*
            (system.branches.A[i][:,n+2].^0.5.-
            sqrt.(system.branches.A0[i][end]))./mmHgToPa;
    end
end
