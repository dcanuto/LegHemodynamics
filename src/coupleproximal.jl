function coupleproximal!(system::LegSystem,n::Int64,h::Float64,uprox::Float64)
    # update right-running invariant
    system.branches.W1root = 2*uprox-system.branches.W2root;
    # update proximal A, Q w/ invariants
    system.branches.A[1][1,n+2] = ((2*system.solverparams.rho/
        system.branches.beta[1][end])^2*(0.125*(system.branches.W1root-
        system.branches.W2root) + system.branches.c0[1][end])^4);
    system.branches.Q[1][1,n+2] = (system.branches.A[1][1,n+2]*0.5*(
        system.branches.W1root + system.branches.W2root));
end
