function interiorbcs!(system::LegSystem,n::Int64,split::Int64,
    Q::Vector{Float64},A::Vector{Float64},children::Vector{Int64})
    # children = system.branches.children[split];
    # println("Interior boundary function assignment:")
    if length(children) == 1
        f = LegModule.fsingle!;
        J = LegModule.Jsingle!;
    elseif length(children) == 2
        f = LegModule.fdouble!;
        J = LegModule.Jdouble!;
    elseif length(children) == 3
        f = LegModule.ftriple!;
        J = LegModule.Jtriple!;
    end
    # println("Interior boundary allocation:")
    W = zeros(length(children)+1);
    beta = zeros(length(children)+1);
    A0 = zeros(length(children)+1);
    c0 = zeros(length(children)+1);
    beta[1] = system.branches.beta[split][end];
    A0[1] = system.branches.A0[split][end];
    c0[1] = system.branches.c0[split][end];
    for j = 1:length(children)
        beta[j+1] = system.branches.beta[children[j]][end];
        A0[j+1] = system.branches.A0[children[j]][end];
        c0[j+1] = system.branches.c0[children[j]][end];
    end
    W[1] = system.branches.W1[split];
    W[2:end] = system.branches.W2[children];
    Qnew = zeros(length(children)+1);
    Anew = zeros(length(children)+1);
    iters = Int64[0];
    # println("Interior Newton solver:")
    LegModule.solvesplits!(iters,Qnew,Anew,children,Q,A,W,beta,A0,c0,
        system.solverparams.rho,f,J,system.solverparams.maxiter,
        system.solverparams.maxval,system.solverparams.epsJ,
        system.solverparams.epsN)
    # println("Interior boundary update save:")
    system.solverparams.totaliter += iters[1];
    system.branches.Q[split][system.solverparams.JL,n+2] = Qnew[1];
    system.branches.A[split][system.solverparams.JL,n+2] = Anew[1];
    for j = 1:length(children)
        system.branches.Q[children[j]][1,n+2] = Qnew[j+1];
        system.branches.A[children[j]][1,n+2] = Anew[j+1];
    end
end
