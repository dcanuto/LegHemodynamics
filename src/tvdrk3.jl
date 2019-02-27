function tvdrk3!(system::LegSystem,times::CVTimer,n::Int64,splits::Vector{Int64},terms::Vector{Int64},uprox::Float64)
    # allocation for invariant/boundary value RK3 loop
    dW = zeros(2);
    W10 = zeros(2,length(system.branches.ID));
    W20 = zeros(2,length(system.branches.ID));
    W11 = zeros(2,length(system.branches.ID));
    W21 = zeros(2,length(system.branches.ID));
    W12 = zeros(2,length(system.branches.ID));
    W22 = zeros(2,length(system.branches.ID));
    lf = zeros(length(system.branches.ID));
    lb = zeros(length(system.branches.ID));
    ret1 = zeros(2);
    ret2 = zeros(2);
    ret3 = zeros(1);
    ret4 = zeros(1);

    # allocation for interior value RK3 loop
    dU = zeros((length(system.solverparams.colsint))*2);
    Q1 = zeros(length(system.solverparams.colsint),length(system.branches.ID));
    A1 = zeros(length(system.solverparams.colsint),length(system.branches.ID));
    Q2 = zeros(length(system.solverparams.colsint),length(system.branches.ID));
    A2 = zeros(length(system.solverparams.colsint),length(system.branches.ID));
    Acat = zeros(length(system.solverparams.colsint)+2);
    Qcat = zeros(length(system.solverparams.colsint)+2);

    # allocation for WENO3 spatial derivative subroutine
    ISA = zeros(3);
    ISQ = zeros(3);
    PA = Vector{Function}(length(system.solverparams.colsint));
    PQ = Vector{Function}(length(system.solverparams.colsint));
    Aint = Vector{Float64}(length(system.solverparams.colsint));
    Qint = Vector{Float64}(length(system.solverparams.colsint));
    wA = Vector{Float64}(3);
    wQ = Vector{Float64}(3);
    c = Vector{Float64}(3);
    a = Vector{Float64}(2);
    p = Vector{Float64}(2);
    x = Vector{Float64}(2);
    Al = Vector{Float64}(2);
    Ar = Vector{Float64}(2);
    Ql = Vector{Float64}(2);
    Qr = Vector{Float64}(2);
    Fl = Vector{Float64}(2);
    Fr = Vector{Float64}(2);
    Hl = Vector{Float64}(2);
    Hr = Vector{Float64}(2);
    Pi = Vector{Float64}(3);
    σ = Vector{Float64}(2);
    eigs = Vector{Float64}(2);
    abso = Vector{Float64}(2);
    J = zeros(2,2);
    Pret = zeros(1);
    s = zeros(1);

    # first RK step
    # println("First RK:")
    tic();
    for i = 1:length(system.branches.ID)
        # println("First step, artery $i, time step $n.")
        # scalings
        # println("Scalings:")
        As = ((system.solverparams.rho.*system.branches.c0[i][end].^2)./system.branches.beta[i][end]).^2;
        Qs = As.*system.branches.c0[i][end];
        xs = system.branches.lengthinmm[i].*LegHemodynamics.mmTom .- system.branches.k[i];
        ts = As.*xs./Qs;
        Acat .= system.branches.A[i][:,n+1]/As;
        Qcat .= system.branches.Q[i][:,n+1]/Qs;
        # println("Initial area: $(system.branches.A[i][:,n+1]./As)")
        # println("Initial flowrate: $(system.branches.Q[i][:,n+1]./Qs)")
        # invariant step
        # println("Invariant calculation:")
        invariants!(ret1,ret2,ret3,ret4,Acat,Qcat,
            system.branches.c0[i][end],system.branches.beta[i][end],system.solverparams.rho,
            system.solverparams.mu,system.solverparams.diffusioncoeff,Qs,As)
        W10[:,i] .= ret1;
        W20[:,i] .= ret2;
        lf[i] = ret3[1];
        lb[i] = ret4[1];
        # println("Forward invariant at time step $n, artery $i: $(W1n[i]) m/s")
        # println("Invariant update:")
        invariantodes!(dW,W10[:,i],W20[:,i],lf[i],lb[i],Qcat,Acat,
            system.solverparams.rho,system.solverparams.mu,system.solverparams.diffusioncoeff,
            system.branches.k[i]./xs,ts,xs,Qs,As);
        # println("Invariant assignment:")
        if i == 1 # aortic root
            system.branches.W1[i] = (W10[end,i] .+ system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2root = (W20[1,i] .+ system.solverparams.h./ts.*dW[2]).*Qs./As;
        elseif isempty(system.branches.children[i]) # terminals
            system.branches.W1end[i] = (W10[end,i] .+ system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] = (W20[1,i] .+ system.solverparams.h./ts.*dW[2]).*Qs./As;
        else # all others
            system.branches.W1[i] = (W10[end,i] .+ system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] = (W20[1,i] .+ system.solverparams.h./ts.*dW[2]).*Qs./As;
        end
        # primitive step
        # println("WENO3:")
        weno3!(dU,system.branches.k[i]./xs,Acat,
            Qcat,system.branches.beta[i][end],
            system.solverparams.rho,system.solverparams.mu,
            system.solverparams.diffusioncoeff,ts,xs,Qs,As,ISA,ISQ,PA,PQ,wA,wQ,
            c,a,p,x,Al,Ar,Ql,Qr,Hl,Hr,Fl,Fr,Pi,σ,eigs,abso,J,Aint,Qint,Pret,s)
        # println("Intermediate step:")
        A1[:,i] .= (system.branches.A[i][2:end-1,n+1]./As .+
            system.solverparams.h./ts.*dU[system.solverparams.acols[1:end-2]]).*As;
        Q1[:,i] .= (system.branches.Q[i][2:end-1,n+1]./Qs .+
            system.solverparams.h./ts.*dU[system.solverparams.qcols[1:end.-2].-2]).*Qs;
    end
    times.tfd += toq();

    # first set of boundary values
    # println("Split solve:")
    tic();
    for i = 1:length(splits)
        # println("Solving interior boundary, parent ID: $(splits[i])")
        children = system.branches.children[splits[i]];
        Q = zeros(length(children)+1);
        A = zeros(length(children)+1);
        Q[1] = system.branches.Q[splits[i]][system.solverparams.JL,n+1];
        A[1] = system.branches.A[splits[i]][system.solverparams.JL,n+1];
        for j = 1:length(children)
            Q[j+1] = system.branches.Q[children[j]][1,n+1];
            A[j+1] = system.branches.A[children[j]][1,n+1];
        end
        LegHemodynamics.interiorbcs!(system,n,splits[i],Q,A,children);
    end
    times.ts += toq();
    # println("0D-1D coupling:")
    tic();
    LegHemodynamics.coupling!(system,n,system.solverparams.h,terms,uprox);
    times.tc += toq();

    # second RK step
    # println("Second RK:")
    tic()
    for i = 1:length(system.branches.ID)
        # println("Second step, artery $i, time step $n.")
        Acat[1] = system.branches.A[i][1,n+2];
        Acat[2:end-1] .= A1[:,i];
        Acat[end] = system.branches.A[i][end,n+2];
        Qcat[1] = system.branches.Q[i][1,n+2];
        Qcat[2:end-1] .= Q1[:,i];
        Qcat[end] = system.branches.Q[i][end,n+2];
        # scalings
        As = ((system.solverparams.rho.*system.branches.c0[i][end].^2)./system.branches.beta[i][end]).^2;
        Qs = As.*system.branches.c0[i][end];
        xs = system.branches.lengthinmm[i].*LegHemodynamics.mmTom .- system.branches.k[i];
        ts = As.*xs./Qs;
        Acat /= As;
        Qcat /= Qs;
        # println("Area after first step: $(Acat./As)")
        # println("Flowrate after first step: $(Qcat./Qs)")
        # invariant step
        invariants!(ret1,ret2,ret3,ret4,Acat,Qcat,
            system.branches.c0[i][end],system.branches.beta[i][end],system.solverparams.rho,
            system.solverparams.mu,system.solverparams.diffusioncoeff,Qs,As)
        W11[:,i] .= ret1;
        W21[:,i] .= ret2;
        lf[i] = ret3[1];
        lb[i] = ret4[1];
        invariantodes!(dW,W11[:,i],W21[:,i],lf[i],lb[i],Qcat,Acat,
            system.solverparams.rho,system.solverparams.mu,system.solverparams.diffusioncoeff,
            system.branches.k[i]./xs,ts,xs,Qs,As);
        if i == 1 # aortic root
            system.branches.W1[i] = (0.75.*W10[end,i] .+ 0.25.*W11[end,i] .+
                0.25.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2root = (0.75.*W20[1,i] .+ 0.25.*W21[1,i] .+
                0.25.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        elseif isempty(system.branches.children[i]) # terminals
            system.branches.W1end[i] = (0.75.*W10[end,i] .+ 0.25.*W11[end,i] .+
                0.25.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] = (0.75.*W20[1,i] .+ 0.25.*W21[1,i] .+
                0.25.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        else # all others
            system.branches.W1[i] = (0.75.*W10[end,i] .+ 0.25.*W11[end,i] .+
                0.25.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] = (0.75.*W20[1,i] .+ 0.25.*W21[1,i] .+
                0.25.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        end
        # primitive step
        weno3!(dU,system.branches.k[i]./xs,Acat,
            Qcat,system.branches.beta[i][end],
            system.solverparams.rho,system.solverparams.mu,
            system.solverparams.diffusioncoeff,ts,xs,Qs,As,ISA,ISQ,PA,PQ,wA,wQ,
            c,a,p,x,Al,Ar,Ql,Qr,Hl,Hr,Fl,Fr,Pi,σ,eigs,abso,J,Aint,Qint,Pret,s)
        A2[:,i] .= (0.75.*system.branches.A[i][2:end-1,n+1]./As .+ 0.25.*A1[:,i]./As .+
            0.25.*system.solverparams.h./ts.*dU[system.solverparams.acols[1:end-2]]).*As;
        Q2[:,i] .= (0.75.*system.branches.Q[i][2:end-1,n+1]./Qs .+ 0.25.*Q1[:,i]./Qs .+
            0.25.*system.solverparams.h./ts.*dU[system.solverparams.qcols[1:end.-2].-2]).*Qs;
    end
    times.tfd += toq();

    # second set of boundary values
    tic();
    for i = 1:length(splits)
        children = system.branches.children[splits[i]];
        Q = zeros(length(children)+1);
        A = zeros(length(children)+1);
        Q[1] = system.branches.Q[splits[i]][system.solverparams.JL,n+2];
        A[1] = system.branches.A[splits[i]][system.solverparams.JL,n+2];
        for j = 1:length(children)
            Q[j+1] = system.branches.Q[children[j]][1,n+2];
            A[j+1] = system.branches.A[children[j]][1,n+2];
        end
        LegHemodynamics.interiorbcs!(system,n,splits[i],Q,A,children);
    end
    times.ts += toq();
    tic()
    LegHemodynamics.coupling!(system,n,0.5*system.solverparams.h,terms,uprox);
    times.tc += toq();

    # third RK step
    # println("Third RK:")
    tic();
    for i = 1:length(system.branches.ID)
        # println("Third step, artery $i, time step $n.")
        Acat[1] = system.branches.A[i][1,n+2];
        Acat[2:end-1] .= A2[:,i];
        Acat[end] = system.branches.A[i][end,n+2];
        Qcat[1] = system.branches.Q[i][1,n+2];
        Qcat[2:end-1] .= Q2[:,i];
        Qcat[end] = system.branches.Q[i][end,n+2];
        # scalings
        As = ((system.solverparams.rho.*system.branches.c0[i][end].^2)./system.branches.beta[i][end]).^2;
        Qs = As.*system.branches.c0[i][end];
        xs = system.branches.lengthinmm[i].*LegHemodynamics.mmTom .- system.branches.k[i];
        ts = As.*xs./Qs;
        Acat /= As;
        Qcat /= Qs;
        # println("Area after second step: $(Acat./As)")
        # println("Flowrate after second step: $(Qcat./Qs)")
        # invariant step
        invariants!(ret1,ret2,ret3,ret4,Acat,Qcat,
            system.branches.c0[i][end],system.branches.beta[i][end],system.solverparams.rho,
            system.solverparams.mu,system.solverparams.diffusioncoeff,Qs,As)
        W12[:,i] .= ret1;
        W22[:,i] .= ret2;
        lf[i] = ret3[1];
        lb[i] = ret4[1];
        invariantodes!(dW,W12[:,i],W22[:,i],lf[i],lb[i],Qcat,Acat,
            system.solverparams.rho,system.solverparams.mu,system.solverparams.diffusioncoeff,
            system.branches.k[i]./xs,ts,xs,Qs,As);
        if i == 1 # aortic root
            system.branches.W1[i] .= (1/3.*W10[end,i] .+ 2/3.*W12[end,i] .+
                2/3.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2root = (1/3.*W20[1,i] .+ 2/3.*W22[1,i] .+
                2/3.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        elseif isempty(system.branches.children[i]) # terminals
            system.branches.W1end[i] .= (1/3.*W10[end,i] .+ 2/3.*W12[end,i] .+
                2/3.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] .= (1/3.*W20[1,i] .+ 2/3.*W22[1,i] .+
                2/3.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        else # all others
            system.branches.W1[i] .= (1/3.*W10[end,i] .+ 2/3.*W12[end,i] .+
                2/3.*system.solverparams.h./ts.*dW[1]).*Qs./As;
            system.branches.W2[i] .= (1/3.*W20[1,i] .+ 2/3.*W22[1,i] .+
                2/3.*system.solverparams.h./ts.*dW[2]).*Qs./As;
        end
        # primitive step
        weno3!(dU,system.branches.k[i]./xs,Acat,
            Qcat,system.branches.beta[i][end],
            system.solverparams.rho,system.solverparams.mu,
            system.solverparams.diffusioncoeff,ts,xs,Qs,As,ISA,ISQ,PA,PQ,wA,wQ,
            c,a,p,x,Al,Ar,Ql,Qr,Hl,Hr,Fl,Fr,Pi,σ,eigs,abso,J,Aint,Qint,Pret,s)
        system.branches.A[i][2:end-1,n+2] .= (1/3.*system.branches.A[i][2:end-1,n+1]./As .+ 2/3.*A2[:,i]./As .+
            2/3.*system.solverparams.h./ts.*dU[system.solverparams.acols[1:end-2]]).*As;
        system.branches.Q[i][2:end-1,n+2] .= (1/3.*system.branches.Q[i][2:end-1,n+1]./Qs .+ 2/3.*Q2[:,i]./Qs .+
            2/3.*system.solverparams.h./ts.*dU[system.solverparams.qcols[1:end.-2].-2]).*Qs;
    end
    times.tfd += toq();

    # third set of boundary values
    # println("Third interior boundaries:")
    tic();
    for i = 1:length(splits)
        # println("Interior boundary allocation:")
        children = system.branches.children[splits[i]];
        Q = zeros(length(children)+1);
        A = zeros(length(children)+1);
        Q[1] = system.branches.Q[splits[i]][system.solverparams.JL,n+2];
        A[1] = system.branches.A[splits[i]][system.solverparams.JL,n+2];
        for j = 1:length(children)
            Q[j+1] = system.branches.Q[children[j]][1,n+2];
            A[j+1] = system.branches.A[children[j]][1,n+2];
        end
        # println("Interior boundary solution:")
        LegHemodynamics.interiorbcs!(system,n,splits[i],Q,A,children);
    end
    times.ts += toq();
    # println("Third coupling:")
    tic();
    LegHemodynamics.coupling!(system,n,system.solverparams.h,terms,uprox);
    times.tc += toq();

end
