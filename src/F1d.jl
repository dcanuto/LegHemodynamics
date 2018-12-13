function F1d(F::Vector{Float64},Q::Float64,A::Float64,beta::Float64,rho::Float64,ts::Float64,
    xs::Float64,Qs::Float64,As::Float64)
    F[1] = ts*Qs./(As*xs)*Q;
    F[2] = Qs*ts/(As*xs)*(Q^2)/A + ts*(As^1.5)./(Qs*xs)*beta/(3*rho)*A^1.5;
    # println("A flux: $(ts*Qs./(As*xs)*Q)")
    # println("Q flux, term 1: $(Qs*ts/(As*xs)*Q^2/A)")
    # println("Q flux, term 2: $(ts*As^1.5./(Qs*xs)*beta/(3*rho)*A^1.5)")
end
