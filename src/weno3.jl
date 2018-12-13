function weno3!(dU::Vector{Float64},k::Float64,A::Vector{Float64},Q::Vector{Float64},
    beta::Float64,rho::Float64,μ::Float64,diff::Float64,ts::Float64,xs::Float64,
    Qs::Float64,As::Float64,ISA::Vector{Float64},ISQ::Vector{Float64},
    PA::Vector{Function},PQ::Vector{Function},wA::Vector{Float64},wQ::Vector{Float64},
    c::Vector{Float64},a::Vector{Float64},p::Vector{Float64},x::Vector{Float64},
    Al::Vector{Float64},Ar::Vector{Float64},Ql::Vector{Float64},Qr::Vector{Float64},
    Hl::Vector{Float64},Hr::Vector{Float64},Fl::Vector{Float64},Fr::Vector{Float64},
    Pi::Vector{Float64},σ::Vector{Float64},eigs::Vector{Float64},abso::Vector{Float64},
    J::Matrix{Float64},Aint::Vector{Float64},Qint::Vector{Float64},Pret::Vector{Float64},
    s::Vector{Float64})

    # parameters
    ϵ = k^2;

    # allocation
    Aint[:] .= A[2:end-1];
    Qint[:] .= Q[2:end-1];

    # polynomial reconstruction loop (interior points only)
    # println("Reconstruction functions:")
    for j = 1:length(Aint)
        if j == 1
            c[1] = ϵ;
            c[2] = 0.25;
            c[3] = 0.75-ϵ;
        elseif j == length(A)
            c[1] = ϵ;
            c[2] = 0.25;
            c[3] = 0.75-ϵ
        else
            c[1] = 0.25;
            c[2] = 0.5;
            c[3] = 0.25;
        end
        smoothinds!(ISA,Aint,c,k,j,σ)
        smoothinds!(ISQ,Qint,c,k,j,σ)
        weights!(wA,c,ϵ,ISA,s)
        weights!(wQ,c,ϵ,ISQ,s)
        PA[j] = reconstruct(wA,c,Aint,j,k,σ)
        PQ[j] = reconstruct(wQ,c,Qint,j,k,σ)
    end

    # rhs loop
    # println("RHS loop:")
    for j = 1:length(Aint)
        if j == 1 # leftmost interior point
            # cell boundaries
            x[1] = 0.;
            x[2] = k;
            # reconstructed values at cell boundaries
            Al[1] = A[1];
            PA[j](Pret,x[1],Pi);
            Al[2] = Pret[1];
            PA[j](Pret,x[2],Pi);
            Ar[1] = Pret[1];
            PA[j+1](Pret,x[2],Pi);
            Ar[2] = Pret[1];
            Ql[1] = Q[1];
            PQ[j](Pret,x[1],Pi);
            Ql[2] = Pret[1];
            PQ[j](Pret,x[2],Pi);
            Qr[1] = Pret[1];
            PQ[j+1](Pret,x[2],Pi);
            Qr[2] = Pret[1];
            # numerical viscosity coefficient based on spectral radius of flux Jacobian
            abseigs1d(abso,eigs,J,Ql[1],Al[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Ql[2],Al[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[1] = maximum(p);
            abseigs1d(abso,eigs,J,Qr[1],Ar[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Qr[2],Ar[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[2] = maximum(p);
        elseif j == length(Aint) # rightmost interior point
            x[1] = k*(length(A)-3.);
            x[2] = k*(length(A)-2.);

            PA[j-1](Pret,x[1],Pi);
            Al[1] = Pret[1];
            PA[j](Pret,x[1],Pi);
            Al[2] = Pret[1];
            PA[j](Pret,x[2],Pi);
            Ar[1] = Pret[1];
            Ar[2] = A[end];
            PQ[j-1](Pret,x[1],Pi);
            Ql[1] = Pret[1];
            PQ[j](Pret,x[1],Pi);
            Ql[2] = Pret[1];
            PQ[j](Pret,x[2],Pi);
            Qr[1] = Pret[1];
            Qr[2] = Q[end];

            abseigs1d(abso,eigs,J,Ql[1],Al[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Ql[2],Al[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[1] = maximum(p);
            abseigs1d(abso,eigs,J,Qr[1],Ar[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Qr[2],Ar[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[2] = maximum(p);
        else # other interior points
            x[1] = k*(j-1);
            x[2] = k*j;

            PA[j-1](Pret,x[1],Pi);
            Al[1] = Pret[1];
            PA[j](Pret,x[1],Pi);
            Al[2] = Pret[1];
            PA[j](Pret,x[2],Pi);
            Ar[1] = Pret[1];
            PA[j+1](Pret,x[2],Pi);
            Ar[2] = Pret[1];
            PQ[j-1](Pret,x[1],Pi);
            Ql[1] = Pret[1];
            PQ[j](Pret,x[1],Pi);
            Ql[2] = Pret[1];
            PQ[j](Pret,x[2],Pi);
            Qr[1] = Pret[1];
            PQ[j+1](Pret,x[2],Pi);
            Qr[2] = Pret[1];

            abseigs1d(abso,eigs,J,Ql[1],Al[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Ql[2],Al[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[1] = maximum(p);
            abseigs1d(abso,eigs,J,Qr[1],Ar[1],beta,rho,ts,xs,Qs,As);
            p[1] = maximum(abso);
            abseigs1d(abso,eigs,J,Qr[2],Ar[2],beta,rho,ts,xs,Qs,As);
            p[2] = maximum(abso);
            a[2] = maximum(p);
        end
        # numerical fluxes (local Lax-Friedrichs)
        F1d(Fl,Ql[1],Al[1],beta,rho,ts,xs,Qs,As)
        F1d(Fr,Ql[2],Al[2],beta,rho,ts,xs,Qs,As)
        Hl[1] = 0.5*(Fl[1] + Fr[1] - a[1]*(Al[2] - Al[1]))
        Hl[2] = 0.5*(Fl[2] + Fr[2] - a[1]*(Ql[2] - Ql[1]))
        F1d(Fr,Qr[2],Ar[2],beta,rho,ts,xs,Qs,As)
        F1d(Fl,Qr[1],Ar[1],beta,rho,ts,xs,Qs,As)
        Hr[1] = 0.5*(Fl[1] + Fr[1] - a[2]*(Ar[2] - Ar[1]))
        Hr[2] = 0.5*(Fl[2] + Fr[2] - a[2]*(Qr[2] - Qr[1]))
        # 3rd-order WENO spatial derivatives (plus source terms)
        dU[j] = -1./k*(Hr[1] - Hl[1]);
        dU[length(Aint)+j] = -1./k*(Hr[2] - Hl[2]) - ts/As*diff*π*μ/rho*Q[j+1]/A[j+1];
    end
    # println("Calculated dU: $dU")
end
