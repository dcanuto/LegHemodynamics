function applycustomics!(system::LegSystem)
    Pstart = rand(Distributions.Normal(70,0.1))*mmHgToPa;
    for i in 1:length(system.branches.ID)
        system.branches.A[i][:,1] = (Pstart/system.branches.beta[i][end] +
            sqrt(system.branches.A0[i][end]))^2;
        system.branches.Q[i][:,1] = 0;
        system.branches.P[i][:,1] = Pstart/mmHgToPa;
    end
    return system
end
