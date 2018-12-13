function coupledistal!(system::LegSystem,n::Int64,h::Float64,terms::Vector{Int64})
    # yout = zeros(4);
    # for i = 1:length(terms)
    #     # yout = zeros(4);
    #     Vs = system.branches.term[terms[i]].V[n+1,1];
    #     vs = system.branches.c0[terms[i]][end];
    #     ts = system.branches.k[terms[i]]/vs;
    #     iters = Int64[0];
    #     CVModule.newtondist!(yout,iters,Vs,vs,ts,system.branches.term[terms[i]].V[n+1,1],
    #         system.branches.term[terms[i]].C[1],system.branches.term[terms[i]].Q[n+1,1],
    #         system.branches.term[terms[i]].V0[1],system.branches.beta[terms[i]][end],
    #         system.branches.c0[terms[i]][end],system.branches.A0[terms[i]][end],
    #         system.branches.W1end[terms[i]],system.solverparams.rho,CVModule.fdist!,
    #         CVModule.Jdist!,system.solverparams.maxiter,system.solverparams.epsJ,
    #         system.solverparams.epsN,system.solverparams.maxval,h);
    #     system.solverparams.totaliter += iters[1];
    #     system.branches.term[terms[i]].V[n+2,1] = yout[1];
    #     system.branches.term[terms[i]].P[n+2,1] = yout[2];
    #     system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[3];
    #     system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[4];
    # end
end
