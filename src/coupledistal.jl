function coupledistal!(system::LegSystem,n::Int64,h::Float64,terms::Vector{Int64})
    # # W1 at next time step
    # for i = 1:length(terms)
    #     ret1 = Float64[0];
    #     CVModule.endinvariants!(ret1,system.branches.Q[terms[i]][system.solverparams.JL,n+1],
    #         system.branches.A[terms[i]][system.solverparams.JL,n+1],system.branches.Q[terms[i]][system.solverparams.JL-1,n+1],
    #         system.branches.A[terms[i]][system.solverparams.JL-1,n+1],system.branches.c0[terms[i]][end],
    #         system.branches.beta[terms[i]][end],system.solverparams.rho,system.solverparams.mu,
    #         system.solverparams.diffusioncoeff,system.branches.k[terms[i]],system.solverparams.h);
    #     system.branches.W1end[terms[i]] = ret1[1];
    # end
    yout = zeros(4);
    for i = 1:length(terms)
        Ws = system.branches.c0[terms[i]][end];
        R0 = system.branches.term[terms[i]].R[1];
        R1 = system.branches.term[terms[i]].R[2];
        C = system.branches.term[terms[i]].C[1];
        beta = system.branches.beta[terms[i]][end];
        A0 = system.branches.A0[terms[i]][end];
        c0 = system.branches.c0[terms[i]][end];
        Aest = (system.branches.P[terms[i]][system.solverparams.JL,n+1]*mmHgToPa/beta + sqrt.(A0))^2;
        W2est = system.branches.W1end[terms[i]] + 8*(c0 - sqrt(0.5*beta/system.solverparams.rho)*Aest^0.25);
        # println("Initial estimate for W2, artery $(terms[i]), time step $n: $W2est")
        Pn = R0*R1*C/h*system.branches.Q[terms[i]][system.solverparams.JL,n+1] -
            R1*C/h*system.branches.P[terms[i]][system.solverparams.JL,n+1]*mmHgToPa -
            (1+R1*C/h)*beta*sqrt(A0);
        iters = Int64[0];
        LegHemodynamics.newtondist!(yout,iters,Ws,W2est,Pn,C,R0,R1,beta,c0,A0,
            system.branches.W1end[terms[i]],system.solverparams.rho,
            LegHemodynamics.fdist!,LegHemodynamics.Jdist!,system.solverparams.maxiter,
            system.solverparams.epsJ,system.solverparams.epsN,
            system.solverparams.maxval,h);
        system.solverparams.totaliter += iters[1];
        system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[1];
        system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[2];
        system.branches.term[terms[i]].P[n+2,1] = yout[3];
        system.branches.term[terms[i]].Q[n+2,1] = yout[4];
    end
end
