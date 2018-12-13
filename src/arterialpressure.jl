function arterialpressure!(system::LegSystem,n::Int64)
    for i = 1:length(system.branches.ID)
        system.branches.P[i][:,n+2] .= system.branches.beta[i][end].*
            (system.branches.A[i][:,n+2].^0.5.-
            sqrt.(system.branches.A0[i][end]))./mmHgToPa;
    end
end
